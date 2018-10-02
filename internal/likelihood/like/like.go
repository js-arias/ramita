// Copyright (c) 2018, J. Salvador Arias <jsalarias@csnat.unt.edu.ar>
// All rights reserved.
// Distributed under BSD2 license that can be found in the LICENSE file.

// Package like implements the l.like command,
// i.e. print the likelihood of a tree.
package like

import (
	"fmt"
	"os"

	"github.com/js-arias/biodv/cmdapp"
	"github.com/js-arias/ramita/likelihood"

	"github.com/pkg/errors"
)

var cmd = &cmdapp.Command{
	UsageLine: "l.like [-o|--optimize] [-t|--tree <treefile>] <dataset>",
	Short:     "print the likelihood of a tree",
	Long: `
Command l.like reads a tree in parenthetical format and prints its
negative log likelihood under a simple poisson model. If the tree
does not have explicit branch lengths, a default branch length of
0.01 will be used.

If the option -o, or --optimize, is used, then it will try to
improve branch lengths.

The tree will be read from the standard input, unless the option
-t or --tree is defined with a tree file.

Options are:

    -o
    --optimize
      Try to optimize the current branch lengths to increase the
      likelihood.

    -t <treefile>
    --tree <treefile>
      If defined, the tree will be read from the indicated file,
      instead of the standard input.

    <dataset>
      The phylogenetic data matrix. It is a required option.
	`,
	Run:           run,
	RegisterFlags: register,
}

func init() {
	cmdapp.Add(cmd)
}

var treefile string
var optimize bool

func register(c *cmdapp.Command) {
	c.Flag.StringVar(&treefile, "tree", "", "")
	c.Flag.StringVar(&treefile, "t", "", "")
	c.Flag.BoolVar(&optimize, "optimize", false, "")
	c.Flag.BoolVar(&optimize, "o", false, "")
}

func run(c *cmdapp.Command, args []string) error {
	if len(args) != 1 {
		return errors.Errorf("%s: expecting a dataset filename", c.Name())
	}

	f, err := os.Open(args[0])
	if err != nil {
		return errors.Wrapf(err, "%s: while opening %s", c.Name(), args[0])
	}
	defer f.Close()

	m, err := likelihood.NewMatrix(f)
	if err != nil {
		return errors.Wrapf(err, "%s: when parsing matrix", c.Name())
	}

	tf := os.Stdin
	if treefile != "" {
		tf, err = os.Open(treefile)
		if err != nil {
			return errors.Wrapf(err, "%s: while opening %s", c.Name(), treefile)
		}
		defer tf.Close()
	}

	tr, err := likelihood.ReadTree(tf, m)
	if err != nil {
		return errors.Wrapf(err, "%s: when parsing tree", c.Name())
	}
	if optimize {
		fmt.Printf("# Origina tree -log Likelihood: %.6f\n", -tr.Like())
		tr.Refine()
	}
	fmt.Printf("# Tree -log Likelihood:\n%.6f\n", -tr.Like())

	return nil
}
