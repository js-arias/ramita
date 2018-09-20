// Copyright (c) 2018, J. Salvador Arias <jsalarias@csnat.unt.edu.ar>
// All rights reserved.
// Distributed under BSD2 license that can be found in the LICENSE file.

package parsimony

import (
	"math/rand"
	"sort"
	"strings"
	"testing"

	"github.com/js-arias/ramita/matrix"
)

// Benchmark for a Wagner tree
func BenchmarkWagner(b *testing.B) {
	r := strings.NewReader(dnaBlob)
	s := matrix.NewScanner(r)
	m, _ := NewMatrix(s)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		m.Wagner()
	}
}

// Benchmark for a Wagner tree without storing
// previous state copies.
func BenchmarkWagnerNoCopy(b *testing.B) {
	r := strings.NewReader(dnaBlob)
	s := matrix.NewScanner(r)
	m, _ := NewMatrix(s)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		m.noCopyWagner()
	}
}

// Wagner returns a new tree,
// build with the Wagner algorithm and
// a random addition sequence.
func (m *Matrix) noCopyWagner() *Tree {
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
		Chars: make([]uint8, len(m.Out.Chars)),
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
		Anc:   root,
		Chars: make([]uint8, len(m.Out.Chars)),
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

	// add the remaning terminals
	for _, i := range ls[2:] {
		tr.noCopyAddTerm(terms[i])
	}

	return tr
}

// AddTerm adds a new terminal to the tree.
func (tr *Tree) noCopyAddTerm(tm *Terminal) {
	na := &Node{
		Chars: make([]uint8, len(tm.Chars)),
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
		increDown(na)
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
	tr.Nodes = append(tr.Nodes, na, nt)
}
