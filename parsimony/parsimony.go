// Copyright (c) 2018, J. Salvador Arias <jsalarias@csnat.unt.edu.ar>
// All rights reserved.
// Distributed under BSD2 license that can be found in the LICENSE file.

// Package parsimony implements
// a simple parsimony search.
package parsimony

import (
	"math/rand"
	"sort"
	"time"

	"github.com/js-arias/ramita/matrix"

	"github.com/pkg/errors"
)

func init() {
	rand.Seed(time.Now().UnixNano())
}

// A Matrix is a parsimony ready dataset.
type Matrix struct {
	Out   *Terminal
	Names map[string]*Terminal
}

// IsValid returns true,
// if the matrix is valid,
// i.e. if the number of characters
// is the same on all terminals.
func (m *Matrix) IsValid() bool {
	var n int
	if m.Out == nil {
		for _, t := range m.Names {
			n = len(t.Chars)
			break
		}
	} else {
		n = len(m.Out.Chars)
	}
	for _, t := range m.Names {
		if n != len(t.Chars) {
			return false
		}
	}
	return true
}

// A Terminal is a terminal taxon
// with phylogenetic (character) data.
type Terminal struct {
	Name  string
	Chars []uint8
}

// NewMatrix returns a new matrix
// from a matrix scanner.
func NewMatrix(s *matrix.Scanner) (*Matrix, error) {
	block := -1
	var ct matrix.DataType // character type of the current block
	var nchars, cblock int // number of chars, total and in current block

	var empty, empBlock []uint8 // slice of unknowns

	var bmap map[string]bool // terminals read on the current block

	m := &Matrix{Names: make(map[string]*Terminal)}

	for s.Scan() {
		tx := s.Taxon()
		if tx.Block != block {
			// A new block
			block = tx.Block
			ct = tx.Type

			// add data to taxons not defined in the block
			for n, t := range m.Names {
				if bmap[n] {
					continue
				}
				t.Chars = append(t.Chars, empBlock...)
			}
			empty = append(empty, empBlock...)
			bmap = make(map[string]bool)

			nchars += cblock
			cblock = len(tx.Chars)
			empBlock = make([]uint8, len(tx.Chars))
			for i := range empBlock {
				empBlock[i] = matrix.Unknown(ct)
			}
		}
		if len(tx.Chars) != cblock {
			return nil, errors.Errorf("parsimony: on block %d: taxon %s with wrong number of chars: %d, want %d", block, tx.Name, len(tx.Chars), cblock)
		}
		if bmap[tx.Name] {
			return nil, errors.Errorf("parsimony: on block %d: taxon %s repeated", block, tx.Name)
		}
		bmap[tx.Name] = true
		t := m.Names[tx.Name]
		if t == nil {
			t = &Terminal{
				Name:  tx.Name,
				Chars: append([]uint8{}, empty...),
			}
			m.Names[t.Name] = t
			if m.Out == nil {
				m.Out = t
			}
		}
		t.Chars = append(t.Chars, tx.Chars...)
	}
	if err := s.Err(); err != nil {
		return nil, errors.Wrap(err, "parsimony")
	}

	// check last block
	for n, t := range m.Names {
		if bmap[n] {
			continue
		}
		t.Chars = append(t.Chars, empBlock...)
	}

	if !m.IsValid() {
		return nil, errors.New("parsimony: bad formatted matrix")
	}
	return m, nil
}

// A Node is a node of a phylogenetic tree.
type Node struct {
	Anc         *Node     // Ancestor
	Left, Right *Node     // Descendants of the node
	Term        *Terminal // A Terminal (in case the node is a terminal)
	Chars       []uint8   // Down-pass assignations
	Cost        int       // Cost at this node
	charsCopy   []uint8   // A copy of the down-pass assignation
	costCopy    int       // A copy if the cost
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

// Wagner returns a new tree,
// build with the Wagner algorithm and
// a random addition sequence.
func (m *Matrix) Wagner() *Tree {
	// randomize terminal order
	terms := make(map[int]*Terminal, len(m.Names)-1)
	var ls []int
	for _, t := range m.Names {
		if t == m.Out {
			continue
		}
		v := rand.Int()
		ls = append(ls, v)
		terms[v] = t
	}
	sort.Ints(ls)

	// Add the firts three terminals
	tr := &Tree{}
	root := &Node{
		Chars:     make([]uint8, len(m.Out.Chars)),
		charsCopy: make([]uint8, len(m.Out.Chars)),
	}
	tr.Root = root
	tr.Nodes = append(tr.Nodes, root)
	out := &Node{
		Anc:   root,
		Term:  m.Out,
		Chars: m.Out.Chars,
	}
	tr.Nodes = append(tr.Nodes, out)
	n0 := &Node{
		Anc:       root,
		Chars:     make([]uint8, len(m.Out.Chars)),
		charsCopy: make([]uint8, len(m.Out.Chars)),
	}
	tr.Nodes = append(tr.Nodes, n0)
	root.Left = out
	root.Right = n0

	tm := terms[ls[0]]
	t0 := &Node{
		Anc:   n0,
		Term:  tm,
		Chars: tm.Chars,
	}
	tr.Nodes = append(tr.Nodes, t0)
	tm = terms[ls[1]]
	t1 := &Node{
		Anc:   n0,
		Term:  tm,
		Chars: tm.Chars,
	}
	tr.Nodes = append(tr.Nodes, t1)
	n0.Left = t0
	n0.Right = t1
	increDown(n0)

	// make the copy of assignations and costs
	for _, n := range tr.Nodes {
		if n.Term != nil {
			continue
		}
		copy(n.charsCopy, n.Chars)
		n.costCopy = n.Cost
	}

	// add the remaning terminals
	for _, i := range ls[2:] {
		tr.addTerm(terms[i])
	}

	return tr
}

// AddTerm adds a new terminal to the tree.
func (tr *Tree) addTerm(tm *Terminal) {
	na := &Node{
		Chars:     make([]uint8, len(tm.Chars)),
		charsCopy: make([]uint8, len(tm.Chars)),
	}
	nt := &Node{
		Anc:   na,
		Term:  tm,
		Chars: tm.Chars,
	}
	na.Left = nt

	var bestPos *Node
	bestCost := tr.Cost() + len(tm.Chars)*2
	for _, d := range tr.Nodes[2:] {
		// Test the position
		a := d.Anc
		na.Anc = a
		na.Right = d
		d.Anc = na
		if a.Left == d {
			a.Left = na
		} else {
			a.Right = na
		}

		cost := increDown(na)
		if cost < bestCost {
			bestCost = cost
			bestPos = d
		}

		// Restore the position
		if a.Left == na {
			a.Left = d
		} else {
			a.Right = d
		}
		d.Anc = a

		// Restore the assignations
		for a != nil {
			copy(a.Chars, a.charsCopy)
			a.Cost = a.costCopy
			a = a.Anc
		}
	}

	// Add the nodes
	a := bestPos.Anc
	na.Anc = a
	na.Right = bestPos
	bestPos.Anc = na
	if a.Left == bestPos {
		a.Left = na
	} else {
		a.Right = na
	}
	increDown(na)

	// Set assignations
	for x := na; x != nil; x = x.Anc {
		copy(x.charsCopy, x.Chars)
		x.costCopy = x.Cost
	}
	tr.Nodes = append(tr.Nodes, na, nt)
}

// IncreDown implements a simple incremental downpass,
// that optimize a node and its ancestors.
// It returns the final cost of the optimization.
func increDown(n *Node) int {
	cost := 0
	for n != nil {
		if n.Term != nil {
			n = n.Anc
			continue
		}
		optimize(n)
		cost = n.Cost
		n = n.Anc
	}
	return cost
}

// Optimize makes an optimization,
// of the current node.
func optimize(n *Node) {
	if n.Term != nil {
		return
	}
	n.Cost = n.Left.Cost + n.Right.Cost
	for i := range n.Chars {
		v := n.Left.Chars[i] & n.Right.Chars[i]
		if v == 0 {
			v = n.Left.Chars[i] | n.Right.Chars[i]
			n.Cost++
		}
		n.Chars[i] = v
	}
}