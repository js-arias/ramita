// Copyright (c) 2018, J. Salvador Arias <jsalarias@csnat.unt.edu.ar>
// All rights reserved.
// Distributed under BSD2 license that can be found in the LICENSE file.

package parsimony

import (
	"bufio"
	"fmt"
	"io"
	"strings"
	"unicode"

	"github.com/js-arias/ramita/matrix"

	"github.com/pkg/errors"
)

// A Node is a node of a phylogenetic tree.
type Node struct {
	Anc         *Node            // Ancestor
	Left, Right *Node            // Descendants of the node
	Term        *matrix.Terminal // A Terminal (in case the node is a terminal)
	Chars       []uint8          // Down-pass assignations
	Cost        int              // Cost at this node
	charsCopy   []uint8          // A copy of the down-pass assignation
	costCopy    int              // A copy if the cost
}

// A Tree is a phylogenetic tree.
type Tree struct {
	Root  *Node   // The root node
	Nodes []*Node // A list of nodes
}

// Cost returns the current cost of the tree.
func (t *Tree) Cost() int {
	return t.Root.Cost
}

// Write writes a tree into a io.Writer.
func (t *Tree) Write(w io.Writer, comma bool) {
	t.Root.write(w, comma)
	fmt.Fprintf(w, ";")
}

// Write write a node into a io.Writer.
func (n *Node) write(w io.Writer, comma bool) {
	if n.Term != nil {
		fmt.Fprintf(w, "%s", n.Term.Name)
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
}

// Laderize moves smaller branches to be left descendants,
// or to be right descendants if right is true.
func (t *Tree) Laderize(right bool) {
	t.Root.laderize(right)
}

// Laderize moves smaller branches to be left descendants,
// or to be right descendants if right is true.
// It returns the number of terminals in the branch.
func (n *Node) laderize(right bool) int {
	if n.Term != nil {
		return 1
	}
	l := n.Left.laderize(right)
	r := n.Right.laderize(right)
	if right {
		if r > l {
			n.Left, n.Right = n.Right, n.Left
		}
		return l + r
	}
	if l > r {
		n.Left, n.Right = n.Right, n.Left
	}
	return l + r
}

// ReadTree reads a tree from a Reader.
func ReadTree(in io.Reader, m *matrix.Matrix) (*Tree, error) {
	r := bufio.NewReader(in)
	for {
		r1, _, err := r.ReadRune()
		if err != nil {
			return nil, errors.Wrapf(err, "parsimony: readtree: unable to read tree")
		}
		if r1 == '(' {
			break
		}
	}
	tr := &Tree{}
	terms := make(map[string]bool)
	root, err := tr.readNode(r, nil, m, terms)
	if err != nil {
		return nil, errors.Wrap(err, "parsimony: readtree")
	}
	tr.Root = root
	return tr, nil
}

// ReadNode reads a node from an reader.
func (tr *Tree) readNode(r *bufio.Reader, anc *Node, m *matrix.Matrix, terms map[string]bool) (*Node, error) {
	n := &Node{Anc: anc}
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
		if r1 == ':' {
			skipBrLen(r)
			continue
		}
		if r1 == ')' {
			break
		}
		if r1 == '(' {
			d, err := tr.readNode(r, n, m, terms)
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
		name, err := readTerm(r)
		if err != nil {
			return nil, err
		}
		tm := m.Names[name]
		if tm == nil {
			return nil, errors.Errorf("terminal %s not in matrix", name)
		}
		if terms[name] {
			return nil, errors.Errorf("terminal %s repeated", name)
		}
		terms[name] = true

		nt := &Node{
			Anc:   n,
			Term:  tm,
			Chars: tm.Chars,
		}
		if n.Left == nil {
			n.Left = nt
		} else if n.Right == nil {
			n.Right = nt
		} else {
			return nil, errors.New("polytomic tree")
		}
	}

	if n.Left == nil || n.Right == nil {
		return nil, errors.New("node without two descendants")
	}
	n.Chars = make([]uint8, len(n.Left.Chars))
	optimize(n)
	n.charsCopy = make([]uint8, len(n.Chars))
	copy(n.charsCopy, n.Chars)
	n.costCopy = n.Cost
	return n, nil
}

// readTerm reads a terminal name on a tree.
func readTerm(r *bufio.Reader) (string, error) {
	var b strings.Builder
	for {
		r1, _, err := r.ReadRune()
		if err != nil {
			return "", err
		}
		if unicode.IsSpace(r1) {
			break
		}
		if r1 == ':' {
			skipBrLen(r)
			break
		}
		if r1 == ',' || r1 == '(' || r1 == ')' {
			r.UnreadRune()
			break
		}
		b.WriteRune(r1)
	}
	return b.String(), nil
}

// skipBrLen skips branch lengths.
func skipBrLen(r *bufio.Reader) {
	for {
		r1, _, err := r.ReadRune()
		if err != nil {
			return
		}
		if unicode.IsSpace(r1) {
			return
		}
		if r1 == ',' || r1 == '(' || r1 == ')' {
			r.UnreadRune()
			return
		}
	}
}
