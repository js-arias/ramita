// Copyright (c) 2018, J. Salvador Arias <jsalarias@csnat.unt.edu.ar>
// All rights reserved.
// Distributed under BSD2 license that can be found in the LICENSE file.

package likelihood

import (
	"strings"
	"testing"
)

var treeLenBlob = `
(Acanthopleura_japonica:0.034103,(((((((Anemonia_sulcata:0.146466,Branchiostoma_floridae:0.073691):0.024796,(Antedon_serrata:0.075537,Balanoglossus_carnosus:0.060488):0.028735):0.034692,((Aphonopelma_sp.:0.056687,Berndtia_purpurea:0.087972):0.032245,Gordius_aquaticus:0.108648):0.023030):0.016813,Brachionus_plicatilis:0.084194):0.018419,(Eisenia_foetida:0.038490,Enchytraeus_sp.:0.038690):0.045899):0.014294,(Barentsia_hildegardae:0.065558,(Chaetonotus_sp.:0.116416,((Dicyema_sp.:0.189882,Gnathostomula_paradoxa:0.189686):0.051786,((Discocelis_tigrina:0.075342,Geocentrophora_sp.:0.085036):0.015107,(Fasciolopsis_bushi:0.097612,Grillotia_erinaceus:0.142001):0.065829):0.028131):0.022083):0.038127):0.018112):0.016258,(Argopecten_irradians:0.017913,Chlamys_islandica:0.011889):0.045857):0.008526);
`

func TestReadTree(t *testing.T) {
	r := strings.NewReader(dnaBlob)
	m, err := NewMatrix(r)
	if err != nil {
		t.Errorf("likelihood: readtree: unexpected error while reading matrix: %v", err)
	}

	tr, err := ReadTree(strings.NewReader(treeLenBlob), m)
	if err != nil {
		t.Errorf("likelihood: readtree: unexpected error while reading tree: %v", err)
	}

	added := make(map[string]bool)
	nt := checkTerminals(t, tr.Root, added)
	if nt != 21 {
		t.Errorf("parsinomy: readtree: tree size %d terminals, want %d", nt, 21)
	}
	for nm := range m.M.Names {
		if !added[nm] {
			t.Errorf("likelihood: readtree: taxon %s not added", nm)
		}
	}
}

func checkTerminals(t *testing.T, n *Node, added map[string]bool) int {
	if n.Term != nil {
		added[n.Term.Name] = true
		return 1
	}
	return checkTerminals(t, n.Left, added) + checkTerminals(t, n.Right, added)
}
