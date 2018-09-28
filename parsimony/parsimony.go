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
)

func init() {
	rand.Seed(time.Now().UnixNano())
}

// Wagner returns a new tree,
// build with the Wagner algorithm and
// a random addition sequence.
func Wagner(m *matrix.Matrix) *Tree {
	// randomize terminal order
	terms := make(map[int]*matrix.Terminal, len(m.Names)-1)
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
func (tr *Tree) addTerm(tm *matrix.Terminal) {
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

		cost, stop := increBound(na, bestCost)
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
			if a == stop {
				break
			}
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

// IncreBound implements a simple incremental downpass,
// and stopped the optimization
// when the cost is greather than a given bound.
// It returns the final cost,
// and the stopping node.
func increBound(n *Node, bound int) (int, *Node) {
	cost := 0
	for n != nil {
		if n.Term != nil {
			n = n.Anc
			continue
		}
		optimize(n)
		cost = n.Cost
		if cost > bound {
			return cost, n
		}
		n = n.Anc
	}
	return cost, nil
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

// Dayoff performs an SPR branch swapping
// on a tree.
func (tr *Tree) Dayoff() {
	// randomize node order
	nodes := make(map[int]*Node, len(tr.Nodes))
	ls := make([]int, 0, len(tr.Nodes))
	for _, n := range tr.Nodes {
		copy(n.charsCopy, n.Chars)
		n.costCopy = n.Cost
		if n == tr.Root {
			continue
		}
		v := rand.Int()
		ls = append(ls, v)
		nodes[v] = n
	}
	sort.Ints(ls)
	for improve := true; improve; {
		improve = tr.swap(nodes, ls)
	}
}

// Swap test a node position among all
// nodes in the indicated node set.
// It returns true if a new position is found.
func (tr *Tree) swap(nodes map[int]*Node, ls []int) bool {
	improved := false
	bestCost := tr.Cost()
	for _, i := range ls {
		imp := false

		// removes the node
		n := nodes[i]
		if n.Anc == tr.Root {
			continue
		}

		a := n.Anc
		sis := a.Left
		if sis == n {
			sis = a.Right
		}
		a.Left = n
		a.Right = nil

		gf := a.Anc
		unc := gf.Left
		if unc == a {
			unc = gf.Right
		}
		gf.Left = unc
		gf.Right = sis
		a.Anc = nil
		sis.Anc = gf

		increDown(gf)
		for x := gf; x != nil; x = x.Anc {
			copy(x.charsCopy, x.Chars)
			x.costCopy = x.Cost
		}

		// test positions of the node
		for _, j := range ls {
			p := nodes[j]
			if p.IsDesc(a) {
				continue
			}
			if p == sis {
				continue
			}
			if p.Anc == tr.Root {
				continue
			}

			pa := p.Anc
			psis := pa.Left
			if psis == p {
				psis = pa.Right
			}
			pa.Left = psis
			pa.Right = a
			p.Anc = a
			a.Right = p
			a.Anc = pa

			cost, stop := increBound(a, bestCost)
			if cost < bestCost {
				// The new position is the best
				// so update backups and break
				for x := a; x != nil; x = x.Anc {
					copy(x.charsCopy, x.Chars)
					x.costCopy = x.Cost
				}
				bestCost = cost
				improved = true
				imp = true
				break
			}

			// restore positions
			p.Anc = pa
			pa.Right = p
			a.Anc = nil
			a.Right = nil

			// Restore assignations
			for x := p; x != nil; x = x.Anc {
				copy(x.Chars, x.charsCopy)
				x.Cost = x.costCopy
				if x == stop {
					break
				}
			}
		}

		if imp {
			// improvement in this node
			continue
		}

		// restore the node
		sis.Anc = a
		a.Right = sis
		a.Anc = gf
		gf.Left = unc
		gf.Right = a
		copy(a.Chars, a.charsCopy)
		a.Cost = a.costCopy
		increDown(gf)
		for x := gf; x != nil; x = x.Anc {
			copy(x.charsCopy, x.Chars)
			x.costCopy = x.Cost
		}
	}
	return improved
}

// IsDesc returns true,
// if the node a,
// is an ancestor of n.
func (n *Node) IsDesc(a *Node) bool {
	for n != nil {
		if n == a {
			return true
		}
		n = n.Anc
	}
	return false
}
