// Copyright (c) 2018, J. Salvador Arias <jsalarias@csnat.unt.edu.ar>
// All rights reserved.
// Distributed under BSD2 license that can be found in the LICENSE file.

// Package matrix reads phylogenetic data
// in the form of a phylogenetic data matrix.
package matrix

import (
	"io"

	"github.com/pkg/errors"
)

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
func NewMatrix(r io.Reader) (*Matrix, error) {
	s := NewScanner(r)

	block := -1
	var ct DataType        // character type of the current block
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
				empBlock[i] = Unknown(ct)
			}
		}
		if len(tx.Chars) != cblock {
			return nil, errors.Errorf("matrix: on block %d: taxon %s with wrong number of chars: %d, want %d", block, tx.Name, len(tx.Chars), cblock)
		}
		if bmap[tx.Name] {
			return nil, errors.Errorf("matrix: on block %d: taxon %s repeated", block, tx.Name)
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
		return nil, errors.Wrap(err, "matrix")
	}

	// check last block
	for n, t := range m.Names {
		if bmap[n] {
			continue
		}
		t.Chars = append(t.Chars, empBlock...)
	}

	if !m.IsValid() {
		return nil, errors.New("matrix: bad formatted matrix")
	}
	return m, nil
}
