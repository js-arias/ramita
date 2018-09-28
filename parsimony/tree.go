// Copyright (c) 2018, J. Salvador Arias <jsalarias@csnat.unt.edu.ar>
// All rights reserved.
// Distributed under BSD2 license that can be found in the LICENSE file.

package parsimony

import (
	"fmt"
	"io"

	"github.com/js-arias/ramita/matrix"
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
