// Copyright (c) 2018, J. Salvador Arias <jsalarias@csnat.unt.edu.ar>
// All rights reserved.
// Distributed under BSD2 license that can be found in the LICENSE file.

package likelihood

import "math"

// A Model is an evolutionary model.
type Model interface {
	// Prob is the probability of change
	// from one state to another,
	// with a given branch length.
	Prob(from, to int, blen float64) float64

	// Freq is the frequency of a given state.
	Freq(s int) float64
}

// Poisson is a generic poisson model.
type Poisson float64

// NewPoisson returns a new poisson model
// with a given number of states.
func NewPoisson(states int) Poisson {
	return Poisson(states)
}

// Prob is the probability of change
// from one state to another,
// with a given branch length.
func (p Poisson) Prob(from, to int, blen float64) float64 {
	s := float64(p)
	if from == to {
		return 1/s + math.Exp(-blen)*(s-1)/s
	}
	return 1/s - math.Exp(-blen)/s
}

// Freq is the frequency of a given state.
// In a poisson model,
// is equal on all states.
func (p Poisson) Freq(s int) float64 {
	return 1 / float64(p)
}

// NewJC returns a simple Poisson model
// for DNA 4 states,
// i.e. the Jukes-Cantor model.
func NewJC() Poisson {
	return Poisson(4)
}
