// Copyright (c) 2018, J. Salvador Arias <jsalarias@csnat.unt.edu.ar>
// All rights reserved.
// Distributed under BSD2 license that can be found in the LICENSE file.

package likelihood

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"
	"unicode"

	"github.com/js-arias/ramita/matrix"

	"github.com/pkg/errors"
)

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
		m := tr.M.Model[i]
		like := float64(0)
		for s, p := range c {
			like += p * m.Freq(s)
		}

		fmt.Printf("%d: %.8f\n", i, like)

		logLike += math.Log(like)
	}
	return logLike
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
		mod := m.Model[i]
		for s := range n.Cond[i] {
			prob := n.Left.condState(mod, i, s) * n.Right.condState(mod, i, s)
			n.Cond[i][s] = prob
		}
	}
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
		n.Cond[i] = make(Conditional, m.states[i])
		if n.Term == nil {
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
		Anc:  anc,
		Cond: make([]Conditional, len(tr.M.Model)),
		Len:  0.01,
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
			Cond: make([]Conditional, len(tr.M.Model)),
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
