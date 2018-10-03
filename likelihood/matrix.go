// Copyright (c) 2018, J. Salvador Arias <jsalarias@csnat.unt.edu.ar>
// All rights reserved.
// Distributed under BSD2 license that can be found in the LICENSE file.

package likelihood

import (
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
		if _, ok := m.mds["mk8"]; !ok {
			m.mds["mk8"] = NewPoisson(8)
		}
		m.model[i] = "mk8"
		m.states[i] = 8
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
