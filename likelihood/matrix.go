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
	M     *matrix.Matrix // the base data matrix
	Model []Model        // the model of each character

	states []int // number of states per character
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
		Model:  make([]Model, len(pm.Kind)),
		states: make([]int, len(pm.Kind)),
	}

	for i, k := range pm.Kind {
		if k == matrix.DNA {
			m.Model[i] = NewJC()
			m.states[i] = 4
			continue
		}
		m.Model[i] = NewPoisson(8)
		m.states[i] = 8
	}
	return m, nil
}
