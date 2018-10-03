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

	// States is the number of states of a model.
	States() int

	// Changes is the number of free change types
	// allowed by the model.
	Changes() int

	// ChangeRate returns the change rate
	// of a given change type.
	ChangeRate(tp int) float64

	// SetChangeRate changes the change rate
	// of a given change type.
	SetChangeRate(tp int, r float64)
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

// States is the number of states of a model.
func (p Poisson) States() int {
	return int(p)
}

// Changes is the number of free change types
// allowed by the model.
// In the case of the Poisson model,
// it is 0,
// all changes are fixed.
func (p Poisson) Changes() int {
	return 1
}

// ChangeRate returns the change rate
// of a given change type.
// In the case of the Poisson model,
// it is 1 / states.
func (p Poisson) ChangeRate(tp int) float64 {
	return 1 / float64(p)
}

// SetChangeRate changes the change rate
// of a given change type.
// As change rates in a poisson model
// are fixed,
// it is ignoder
func (p Poisson) SetChangeRate(tp int, r float64) {}

// NewJC returns a simple Poisson model
// for DNA 4 states,
// i.e. the Jukes-Cantor model.
func NewJC() Poisson {
	return Poisson(4)
}
