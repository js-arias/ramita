// Copyright (c) 2018, J. Salvador Arias <jsalarias@csnat.unt.edu.ar>
// All rights reserved.
// Distributed under BSD2 license that can be found in the LICENSE file.

package matrix

import (
	"strings"
	"testing"
)

func TestNewMatrix(t *testing.T) {
	r := strings.NewReader(dnaBlob)
	m, err := NewMatrix(r)
	if err != nil {
		t.Errorf("parsinomy: newmatrix: unexpected error while reading matrix: %v", err)
	}
	if len(m.Names) != 21 {
		t.Errorf("parsimony: newmatrix: taxons in the matrix: %d, want %d", len(m.Names), 21)
	}
	term := m.Names["Dicyema_sp."]
	if term == nil {
		t.Errorf("parsimony: newmatrix: taxon %s not in the matrix", "Dicyema_sp.")
	}
	if len(term.Chars) != 2555 {
		t.Errorf("parsimony: newmatrix: %d characters, want %d", len(term.Chars), 2555)
	}

	r = strings.NewReader(dnaBlob + "\n" + morphoBlob)
	m, err = NewMatrix(r)
	if err != nil {
		t.Errorf("parsinomy: newmatrix: unexpected error while reading matrix: %v", err)
	}
	if len(m.Names) != 42 {
		t.Errorf("parsimony: newmatrix: taxons in the matrix: %d, want %d", len(m.Names), 42)
	}
	term = m.Names["Astyanax_abramis"]
	if term == nil {
		t.Errorf("parsimony: newmatrix: taxon %s not in the matrix", "Dicyema_sp.")
	}
	if len(term.Chars) != 2922 {
		t.Errorf("parsimony: newmatrix: %d characters, want %d", len(term.Chars), 2922)
	}
}
