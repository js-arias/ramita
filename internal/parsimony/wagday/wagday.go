// Copyright (c) 2018, J. Salvador Arias <jsalarias@csnat.unt.edu.ar>
// All rights reserved.
// Distributed under BSD2 license that can be found in the LICENSE file.

// Package wagday implements the p.wagner command,
// i.e. make a Wagner-Dayoff tree with parsimony.
package parsimony

import (
	"fmt"
	"os"

	"github.com/js-arias/biodv/cmdapp"
	"github.com/js-arias/ramita/matrix"
	"github.com/js-arias/ramita/parsimony"

	"github.com/pkg/errors"
)

var cmd = &cmdapp.Command{
	UsageLine: "p.wagday [-c|--comma] [<dataset>]",
	Short:     "make a Wagner-Dayoff tree with parsimony",
	Long: `
Command p.wagday makes a tree with parsimony using a random addition
sequence. The resulting tree will be printed in the standard output.

By default, the tree will be printed with sister groups separed by
spaces (tnt format). If the option -c or --comma is set, then sister
groups will be separated by commas (,) as in phylip.

Options are:

    -c
    --comma
      If set, sister groups will be separated by commas.

    <dataset>
      The phylogenetic data matrix. If not given explicitly, it will
      be read from the standard input.
	`,
	Run:           run,
	RegisterFlags: register,
}

func init() {
	cmdapp.Add(cmd)
}

var comma bool

func register(c *cmdapp.Command) {
	c.Flag.BoolVar(&comma, "comma", false, "")
	c.Flag.BoolVar(&comma, "c", false, "")
}

func run(c *cmdapp.Command, args []string) error {
	if len(args) > 1 {
		return errors.Errorf("%s: too many arguments", c.Name())
	}

	f := os.Stdin
	if len(args) == 1 {
		var err error
		f, err = os.Open(args[0])
		if err != nil {
			return errors.Wrapf(err, "%s: while opening %s", c.Name(), args[0])
		}
		defer f.Close()
	}

	m, err := parsimony.NewMatrix(matrix.NewScanner(f))
	if err != nil {
		return errors.Wrapf(err, "%s: when parsing matrix", c.Name())
	}

	tr := m.Wagner()
	fmt.Printf("# Wagner Length: %d\n", tr.Cost())
	tr.Dayoff()
	tr.Laderize(false)
	fmt.Printf("# Final Length: %d\n", tr.Cost())
	tr.Write(os.Stdout, comma)
	fmt.Printf("\n")
	return nil
}
