// Copyright (c) 2018, J. Salvador Arias <jsalarias@csnat.unt.edu.ar>
// All rights reserved.
// Distributed under BSD2 license that can be found in the LICENSE file.

package likelihood

import (
	"fmt"
	"io"

	"github.com/js-arias/ramita/matrix"

	"github.com/pkg/errors"
)

// A Matrix is a phylogenetic matrix
// for a likelihood analysis.
type Matrix struct {
	M *matrix.Matrix // the base data matrix

	model  []string         // the model of each character
	mds    map[string]Model // list of models assigned to the matrix
	states []int            // number of states per character
}

// NewMatrix returns a new matrix
// from a reader.
func NewMatrix(r io.Reader) (*Matrix, error) {
	pm, err := matrix.NewMatrix(r)
	if err != nil {
		return nil, errors.Wrap(err, "likelihood")
	}

	m := &Matrix{
		M:      pm,
		model:  make([]string, len(pm.Kind)),
		mds:    make(map[string]Model),
		states: make([]int, len(pm.Kind)),
	}

	for i, k := range pm.Kind {
		if k == matrix.DNA {
			if _, ok := m.mds["jc"]; !ok {
				m.mds["jc"] = NewJC()
			}
			m.model[i] = "jc"
			m.states[i] = 4
			continue
		}
		var states uint8
		for _, tx := range m.M.Names {
			if tx.Chars[i] == 255 {
				continue
			}
			states |= tx.Chars[i]
		}
		max := 2
		for b := uint8(7); b > 0; b-- {
			if states&(1<<b) != 0 {
				max = int(b) + 1
				break
			}
		}
		nm := fmt.Sprintf("mk%d", max)
		if _, ok := m.mds[nm]; !ok {
			m.mds[nm] = NewPoisson(max)
		}
		m.model[i] = nm
		m.states[i] = max
	}
	return m, nil
}

// Model returns the model used to estimate
// the maximum likelihood of a character.
func (m *Matrix) Model(char int) Model {
	nm := m.model[char]
	return m.mds[nm]
}

// Chars returns the number of characters
// in the datamatrix.
func (m *Matrix) Chars() int {
	return len(m.model)
}

// States returns the number of states of a character.
func (m *Matrix) States(char int) int {
	return m.states[char]
}
