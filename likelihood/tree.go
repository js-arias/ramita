// Copyright (c) 2018, J. Salvador Arias <jsalarias@csnat.unt.edu.ar>
// All rights reserved.
// Distributed under BSD2 license that can be found in the LICENSE file.

package likelihood

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"math/rand"
	"sort"
	"strconv"
	"strings"
	"time"
	"unicode"

	"github.com/js-arias/ramita/matrix"

	"github.com/pkg/errors"
)

func init() {
	rand.Seed(time.Now().UnixNano())
}

// A Conditional is the conditional likelihood
// of a set of states.
type Conditional []float64

// A Node is a node of a phylogenetic tree.
type Node struct {
	Anc         *Node            // Ancestor
	Left, Right *Node            // Descendants of the node
	Term        *matrix.Terminal // A Terminal (in case the node is a terminal)
	Cond        []Conditional    // Conditional likelihood of each character
	Len         float64          // Length of the current branch

	// backups
	condCopy []Conditional
}

// A Tree is a phylogenetic tree.
type Tree struct {
	Root  *Node
	Nodes []*Node
	M     *Matrix
}

// Like returns the log likelihood of the tree.
func (tr *Tree) Like() float64 {
	logLike := float64(0)
	for i, c := range tr.Root.Cond {
		m := tr.M.Model(i)
		like := float64(0)
		for s, p := range c {
			like += p * m.Freq(s)
		}
		logLike += math.Log(like)
	}
	return logLike
}

// Write writes a tree into a io.Writer.
func (t *Tree) Write(w io.Writer, comma bool) {
	t.Root.write(w, comma)
	fmt.Fprintf(w, ";")
}

// Write write a node into a io.Writer.
func (n *Node) write(w io.Writer, comma bool) {
	if n.Term != nil {
		fmt.Fprintf(w, "%s:%.6f", n.Term.Name, n.Len)
		return
	}
	fmt.Fprintf(w, "(")
	n.Left.write(w, comma)
	if comma {
		fmt.Fprintf(w, ",")
	} else {
		fmt.Fprintf(w, " ")
	}
	n.Right.write(w, comma)
	fmt.Fprintf(w, ")")
	if n.Anc != nil {
		fmt.Fprintf(w, ":%.6f", n.Len)
	}
}

// CondState calculates the conditional
// of state s on a node for the c character.
func (n *Node) condState(m Model, c, s int) float64 {
	probX := float64(0)
	for x, l := range n.Cond[c] {
		probX += m.Prob(s, x, n.Len) * l
	}
	return probX
}

// Optimeze makes an optimization,
// of the current node.
func (n *Node) optimize(m *Matrix) {
	if n.Term != nil {
		return
	}
	for i := range n.Cond {
		mod := m.Model(i)
		for s := range n.Cond[i] {
			prob := n.Left.condState(mod, i, s) * n.Right.condState(mod, i, s)
			n.Cond[i][s] = prob
		}
	}
}

// FullOpt optimize a node
// and all of its descendants.
func (n *Node) fullOpt(m *Matrix, id string) {
	if n.Term != nil {
		return
	}
	n.Left.fullOpt(m, id)
	n.Right.fullOpt(m, id)

	for i := range n.Cond {
		if m.model[i] != id {
			continue
		}
		md := m.Model(i)
		for s := range n.Cond[i] {
			prob := n.Left.condState(md, i, s) * n.Right.condState(md, i, s)
			n.Cond[i][s] = prob
		}
	}
}

// IncreDown implements a simple incremental downpass,
// that optimize a node and its ancestors.
func increDown(n *Node, m *Matrix) {
	for n != nil {
		if n.Term != nil {
			n = n.Anc
			continue
		}
		n.optimize(m)
		n = n.Anc
	}
}

// Estimate perfomrs a simple estimation
// of the model parameters
// under the current branch lengths.
func (tr *Tree) Estimate() {
	// get the model list
	models := make(map[string]bool)
	for _, id := range tr.M.model {
		models[id] = true
	}

	like := tr.Like()
	for {
		for id := range models {
			tr.estimate(id, 0.1)
		}
		l := tr.Like()
		if math.Abs(like-l) < 0.001 {
			break
		}
		like = l
	}
}

// Estimate estimates change parameters
// in a recursive fashion.
func (tr *Tree) estimate(id string, step float64) {
	if step < 0.001 {
		return
	}
	like := tr.Like()
	md := tr.M.mds[id]

	for tp := 0; tp < md.Changes(); tp++ {
		// move rate up
		ref := true
		up := false
		best := md.ChangeRate(tp)
		for ref {
			ref = false
			b := best + step
			if b >= 1 {
				break
			}
			md.SetChangeRate(tp, b)
			tr.Root.fullOpt(tr.M, id)
			l := tr.Like()
			if l > like {
				like = l
				best = b
				ref = true
				up = true
				continue
			}
		}

		md.SetChangeRate(tp, best)
		tr.Root.fullOpt(tr.M, id)
		if up {
			tr.estimate(id, step/10)
			continue
		}

		// move rate down
		ref = true
		for ref {
			ref = false
			b := best - step
			if b <= 0 {
				break
			}
			md.SetChangeRate(tp, b)
			tr.Root.fullOpt(tr.M, id)
			l := tr.Like()
			if l > like {
				like = l
				best = b
				ref = true
				continue
			}
		}

		md.SetChangeRate(tp, best)
		tr.Root.fullOpt(tr.M, id)
		tr.estimate(id, step/10)
	}
}

// Refine permforms a simple
// branch length refinement of the tree.
func (tr *Tree) Refine() {
	// randomize node order
	nodes := make(map[int]*Node, len(tr.Nodes))
	ls := make([]int, 0, len(tr.Nodes))
	for _, n := range tr.Nodes {
		if n == tr.Root {
			continue
		}
		v := rand.Int()
		ls = append(ls, v)
		nodes[v] = n
	}
	sort.Ints(ls)
	like := tr.Like()
	for {
		for _, i := range ls {
			n := nodes[i]
			if n == tr.Root {
				continue
			}
			tr.refine(n, 0.1)
		}
		tr.Estimate()
		l := tr.Like()
		if math.Abs(like-l) < 0.001 {
			break
		}
		like = l
	}
}

// Refine refines a branch length
// in a recursive fashion.
func (tr *Tree) refine(n *Node, step float64) {
	if step < 0.001 {
		return
	}
	like := tr.Like()
	best := n.Len

	// move branch length value up
	ref := true
	up := false
	for ref {
		ref = false
		b := best + step
		if b > 100 {
			break
		}
		n.Len = b
		increDown(n.Anc, tr.M)
		l := tr.Like()
		if l > like {
			like = l
			best = b
			ref = true
			up = true
			continue
		}
	}

	n.Len = best
	increDown(n.Anc, tr.M)
	if up {
		tr.refine(n, step/10)
		return
	}

	// move branch length value down
	ref = true
	for ref {
		ref = false
		b := best - step
		if b < 0.0001 {
			break
		}
		n.Len = b
		increDown(n.Anc, tr.M)
		l := tr.Like()
		if l > like {
			like = l
			best = b
			ref = true
			continue
		}
	}

	n.Len = best
	increDown(n.Anc, tr.M)
	tr.refine(n, step/10)
}

// ReadTree reads a tree from a Reader.
func ReadTree(in io.Reader, m *Matrix) (*Tree, error) {
	r := bufio.NewReader(in)
	for {
		r1, _, err := r.ReadRune()
		if err != nil {
			return nil, errors.Wrapf(err, "likelihood: readtree: unable to read tree")
		}
		if r1 == '(' {
			break
		}
	}
	tr := &Tree{M: m}
	terms := make(map[string]bool)
	root, err := tr.readNode(r, nil, terms)
	if err != nil {
		return nil, errors.Wrap(err, "likelihood: readtree")
	}
	root.Len = 0
	tr.Root = root
	return tr, nil
}

func (n *Node) initializeConditionals(m *Matrix) {
	for i := range n.Cond {
		md := m.Model(i)
		n.Cond[i] = make(Conditional, md.States())
		if n.Term == nil {
			n.condCopy[i] = make(Conditional, md.States())
			continue
		}
		tm := n.Term
		for b := 0; b < m.states[i]; b++ {
			if tm.Chars[i]&(1<<uint8(b)) != 0 {
				n.Cond[i][b] = 1
			}
		}
	}
}

// ReadNode reads a node from an reader.
func (tr *Tree) readNode(r *bufio.Reader, anc *Node, terms map[string]bool) (*Node, error) {
	n := &Node{
		Anc:      anc,
		Cond:     make([]Conditional, tr.M.Chars()),
		Len:      0.01,
		condCopy: make([]Conditional, tr.M.Chars()),
	}
	n.initializeConditionals(tr.M)
	tr.Nodes = append(tr.Nodes, n)

	for {
		r1, _, err := r.ReadRune()
		if err != nil {
			return nil, err
		}
		if unicode.IsSpace(r1) {
			continue
		}
		if r1 == ',' {
			continue
		}
		if r1 == ')' {
			break
		}
		if r1 == '(' {
			d, err := tr.readNode(r, n, terms)
			if err != nil {
				return nil, err
			}
			if n.Left == nil {
				n.Left = d
			} else if n.Right == nil {
				n.Right = d
			} else {
				return nil, errors.New("polytomic tree")
			}
			continue
		}

		// a terminal
		r.UnreadRune()
		name, l, err := readTerm(r)
		if err != nil {
			return nil, err
		}
		tm := tr.M.M.Names[name]
		if tm == nil {
			return nil, errors.Errorf("terminal %s not in matrix", name)
		}
		if terms[name] {
			return nil, errors.Errorf("terminal %s repeated", name)
		}
		terms[name] = true

		nt := &Node{
			Anc:  n,
			Term: tm,
			Len:  l,
			Cond: make([]Conditional, tr.M.Chars()),
		}
		nt.initializeConditionals(tr.M)
		if n.Left == nil {
			n.Left = nt
		} else if n.Right == nil {
			n.Right = nt
		} else {
			return nil, errors.New("polytomic tree")
		}
		tr.Nodes = append(tr.Nodes, nt)
	}
	if n.Left == nil || n.Right == nil {
		return nil, errors.New("node without two descendants")
	}
	n.optimize(tr.M)
	copy(n.condCopy, n.Cond)

	if anc != nil {
		r1, _, err := r.ReadRune()
		if err != nil {
			return nil, err
		}
		if r1 != ':' {
			r.UnreadRune()
		} else {
			l, err := readBrLen(r)
			if err != nil {
				return nil, errors.Wrap(err, "bad branch length")
			}
			n.Len = l
		}
	}
	return n, nil
}

// readTerm reads a terminal name on a tree.
func readTerm(r *bufio.Reader) (string, float64, error) {
	var b strings.Builder
	l := 0.01
	for {
		r1, _, err := r.ReadRune()
		if err != nil {
			return "", 0, err
		}
		if unicode.IsSpace(r1) {
			break
		}
		if r1 == ':' {
			l, err = readBrLen(r)
			if err != nil {
				return "", 0, errors.Wrapf(err, "on terminal %s: bad branch length")
			}
			break
		}
		if r1 == ',' || r1 == '(' || r1 == ')' {
			r.UnreadRune()
			break
		}
		b.WriteRune(r1)
	}
	return b.String(), l, nil
}

// readBrLen skips branch lengths.
func readBrLen(r *bufio.Reader) (float64, error) {
	var b strings.Builder
	for {
		r1, _, err := r.ReadRune()
		if err != nil {
			return 0, err
		}
		if unicode.IsSpace(r1) {
			break
		}
		if r1 == ',' || r1 == '(' || r1 == ')' {
			r.UnreadRune()
			break
		}
		b.WriteRune(r1)
	}
	return strconv.ParseFloat(b.String(), 64)
}
