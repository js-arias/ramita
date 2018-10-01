// Copyright (c) 2018, J. Salvador Arias <jsalarias@csnat.unt.edu.ar>
// All rights reserved.
// Distributed under BSD2 license that can be found in the LICENSE file.

package parsimony

import (
	"strings"
	"testing"

	"github.com/js-arias/ramita/matrix"
)

var treeBlob = `
(Acanthopleura_japonica ((Chlamys_islandica Argopecten_irradians) (Barentsia_hildegardae ((Enchytraeus_sp. Eisenia_foetida) (((Antedon_serrata Balanoglossus_carnosus) (Anemonia_sulcata Branchiostoma_floridae)) ((Gordius_aquaticus (Berndtia_purpurea Aphonopelma_sp.)) (Brachionus_plicatilis ((Chaetonotus_sp. (Dicyema_sp. Gnathostomula_paradoxa)) ((Discocelis_tigrina Geocentrophora_sp.) (Grillotia_erinaceus Fasciolopsis_bushi))))))))));
`

func TestReadTree(t *testing.T) {
	r := strings.NewReader(dnaBlob)
	m, err := matrix.NewMatrix(r)
	if err != nil {
		t.Errorf("parsinomy: readtree: unexpected error while reading matrix: %v", err)
	}

	tr, err := ReadTree(strings.NewReader(treeBlob), m)
	if err != nil {
		t.Errorf("parsinomy: readtree: unexpected error while reading tree: %v", err)
	}
	if tr.Cost() != 3822 {
		t.Errorf("parsimony: readtree: tree length %d, want %d", tr.Cost(), 3822)
	}

	added := make(map[string]bool)
	nt := checkTerminals(t, tr.Root, added)
	if nt != 21 {
		t.Errorf("parsinomy: readtree: tree size %d terminals, want %d", nt, 21)
	}
	for nm := range m.Names {
		if !added[nm] {
			t.Errorf("parsinomy: readtree: taxon %s not added", nm)
		}
	}
}
