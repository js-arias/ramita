// Copyright (c) 2018, J. Salvador Arias <jsalarias@csnat.unt.edu.ar>
// All rights reserved.
// Distributed under BSD2 license that can be found in the LICENSE file.

// Package matrix reads phylogenetic data
// in the form of a phylogenetic data matrix.
package matrix

import (
	"bufio"
	"io"
	"strings"
	"unicode"

	"github.com/pkg/errors"
)

// A Taxon is a taxon name with phylogenetic (character) data.
type Taxon struct {
	Name  string
	Type  DataType
	Chars []uint8
}

// A DataType is the kind of the read phylogenetic data.
type DataType uint

// DataType kinds.
const (
	Morphology DataType = iota // Morphological data, valid states: 0-7
	DNA                        // DNA data, valid states: A, C, G, T (and IUPAC polymorphics)
)

// String returns the name of a data type.
func (d DataType) String() string {
	switch d {
	case Morphology:
		return "morphology"
	case DNA:
		return "dna"
	}
	return "unknown"
}

// A Scanner reads phylogenetic character data from a reader.
type Scanner struct {
	r     *bufio.Reader
	kind  DataType
	block int
	taxon *Taxon
	err   error
}

// NewScanner returns a scanner that reads from r.
func NewScanner(r io.Reader) *Scanner {
	s := &Scanner{r: bufio.NewReader(r)}
	for {
		r1 := peekRune(s.r)
		if r1 == 0 {
			_, _, err := s.r.ReadRune()
			if err == nil {
				err = errors.New("unexpected error")
			}
			s.err = errors.Wrap(err, "while starting scanner")
			return s
		}
		if r1 == '#' {
			s.r.ReadString('\n')
			continue
		}
		if r1 == '>' {
			s.r.ReadRune()
			kind, err := readDataType(s.r)
			if err != nil {
				s.err = errors.Wrap(err, "while starting scanner")
				return s
			}
			s.kind = kind
			s.block = 1
			break
		}
	}
	return s
}

// Scan moves the scanner to the next taxon.
// If there are no more taxons,
// or an error happens while preparing it,
// it will return false.
// A call to Err should be made to discriminate
// among these options.
//
// Every call to Taxon,
// even the first one,
// must be preceded by a Scan call.
func (s *Scanner) Scan() bool {
	if s.err != nil {
		return false
	}
	for {
		r1 := peekRune(s.r)
		if r1 == 0 {
			_, _, err := s.r.ReadRune()
			if err == nil {
				err = errors.Errorf("block %d: unexpected error", s.block)
			}
			s.err = err
			return false
		}
		if r1 == '#' {
			s.r.ReadString('\n')
			continue
		}
		if r1 == '>' {
			s.r.ReadRune()
			kind, err := readDataType(s.r)
			if err != nil {
				s.err = errors.Wrapf(err, "expecting block: %d", s.block+1)
				return false
			}
			s.kind = kind
			s.block++
			continue
		}
		ln, err := s.r.ReadString('\n')
		if err != nil {
			s.err = err
			return false
		}
		entry := strings.Fields(ln)
		name := entry[0]
		if len(entry) == 1 {
			s.err = errors.Errorf("block %d: taxon %s: no characers", s.block, name)
			return false
		}
		var data []uint8
		dr := bufio.NewReader(strings.NewReader(strings.Join(entry[1:], "")))
		for {
			c, err := readStates(dr, s.kind)
			if err == io.EOF {
				break
			}
			if err != nil {
				s.err = errors.Wrapf(err, "block %d: taxon %s", s.block, name)
				return false
			}
			data = append(data, c)
		}
		s.taxon = &Taxon{Name: name, Type: s.kind, Chars: data}
		return true
	}
}

// Taxon returns the last read Taxon.
//
// Every call to Taxon,
// even the first one,
// must be preceded by a Scan call.
func (s *Scanner) Taxon() *Taxon {
	if s.taxon == nil {
		panic("matrix: scan should be called before taxon")
	}
	tax := s.taxon
	s.taxon = nil
	return tax
}

// Err returns the last error
// found during iteration.
func (s *Scanner) Err() error {
	if s.err == io.EOF {
		return nil
	}
	return s.err
}

func readDataType(r *bufio.Reader) (DataType, error) {
	if err := skipSpaces(r); err != nil {
		return 0, err
	}

	var last rune
	var b strings.Builder
	for {
		r1, _, err := r.ReadRune()
		if err != nil {
			return 0, err
		}
		if unicode.IsSpace(r1) {
			last = r1
			break
		}
		b.WriteRune(r1)
	}
	if last != '\n' {
		r.ReadString('\n')
	}
	tp := strings.ToLower(b.String())
	if tp == "dna" {
		return DNA, nil
	}
	if strings.HasPrefix(tp, "morpho") {
		return Morphology, nil
	}
	return 0, errors.Errorf("unknown data type: %s", tp)
}

func readStates(r *bufio.Reader, kind DataType) (uint8, error) {
	if kind == DNA {
		r1, _, err := r.ReadRune()
		if err != nil {
			return 0, err
		}
		switch r1 {
		case 'A', 'a':
			return 1, nil
		case 'C', 'c':
			return 2, nil
		case 'G', 'g':
			return 4, nil
		case 'T', 't', 'U', 'u':
			return 8, nil
		case 'Y', 'y':
			return 2 | 8, nil // C or T
		case 'R', 'r':
			return 1 | 4, nil // A or G
		case 'W', 'w':
			return 1 | 8, nil // A or T
		case 'S', 's':
			return 2 | 4, nil // C or G
		case 'K', 'k':
			return 4 | 8, nil // G or T
		case 'M', 'm':
			return 1 | 2, nil // A or C
		case 'B', 'b':
			return 2 | 4 | 8, nil // not A
		case 'D', 'd':
			return 1 | 4 | 8, nil // not C
		case 'H', 'h':
			return 1 | 2 | 8, nil // not G
		case 'V', 'v':
			return 1 | 2 | 4, nil // not T
		case 'X', 'x', 'N', 'n', '?', '-', 'O', 'o':
			return 15, nil
		}
		return 0, errors.Errorf("unknown symbol %q", r1)
	}
	if kind != Morphology {
		return 0, io.EOF
	}
	for {
		r1, _, err := r.ReadRune()
		if err != nil {
			return 0, err
		}
		if r1 == '?' || r1 == '-' {
			return 255, nil
		}
		if r1 == '[' || r1 == '(' {
			var c uint8
			for {
				r1, _, err := r.ReadRune()
				if err != nil {
					return 0, errors.Wrap(err, "while reading polymorh")
				}
				if r1 == ')' || r1 == ']' {
					return c, nil
				}
				if !unicode.IsDigit(r1) {
					return 0, errors.Errorf("while reading polymorph: unknown symbol %q", r1)
				}
				c |= 1 << (uint8(r1) - uint8('0'))
			}
		}
		if !unicode.IsDigit(r1) {
			return 0, errors.Errorf("unknown symbol %q", r1)
		}
		return 1 << (uint8(r1) - uint8('0')), nil
	}
}

func peekRune(r *bufio.Reader) rune {
	if err := skipSpaces(r); err != nil {
		return 0
	}
	r1, _, _ := r.ReadRune()
	r.UnreadRune()
	return r1
}

func skipSpaces(r *bufio.Reader) error {
	for {
		r1, _, err := r.ReadRune()
		if err != nil {
			return err
		}
		if !unicode.IsSpace(r1) {
			r.UnreadRune()
			return nil
		}
	}
}
