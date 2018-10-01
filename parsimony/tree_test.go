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
var treeLenBlob = `
(Acanthopleura_japonica:0.123 ((Chlamys_islandica:0.011 Argopecten_irradians:0.012):0.17 (Barentsia_hildegardae:0.08 ((Enchytraeus_sp.:0.14 Eisenia_foetida:0.3246):0.43546 (((Antedon_serrata:0.001 Balanoglossus_carnosus:0.123):0.131 (Anemonia_sulcata:0.324 Branchiostoma_floridae:0.23445):0.234):0.32243 ((Gordius_aquaticus:0.00234 (Berndtia_purpurea:0.456 Aphonopelma_sp.:0.2345234):0.008):0.045 (Brachionus_plicatilis:0.0014 ((Chaetonotus_sp.:0.0456 (Dicyema_sp.:0.0345 Gnathostomula_paradoxa:0.0067):0.0056):0.14 ((Discocelis_tigrina:0.0354 Geocentrophora_sp.:0.045):0.23 (Grillotia_erinaceus:0.785 Fasciolopsis_bushi:0.0123):0.4365):0.0012):0.4346):0.04567):0.00015):0.0014):0.0002):0.882):0.345);
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

	tr, err = ReadTree(strings.NewReader(treeLenBlob), m)
	if err != nil {
		t.Errorf("parsinomy: readtree: unexpected error while reading tree: %v", err)
	}
	added = make(map[string]bool)
	nt = checkTerminals(t, tr.Root, added)
	if nt != 21 {
		t.Errorf("parsinomy: readtree: tree size %d terminals, want %d", nt, 21)
	}
	for nm := range m.Names {
		if !added[nm] {
			t.Errorf("parsinomy: readtree: taxon %s not added", nm)
		}
	}
}
