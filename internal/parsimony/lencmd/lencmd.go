// Copyright (c) 2018, J. Salvador Arias <jsalarias@csnat.unt.edu.ar>
// All rights reserved.
// Distributed under BSD2 license that can be found in the LICENSE file.

// Package lencmd implements the p.len command,
// i.e. print the length of a tree.
package lencmd

import (
	"fmt"
	"os"

	"github.com/js-arias/biodv/cmdapp"
	"github.com/js-arias/ramita/matrix"
	"github.com/js-arias/ramita/parsimony"

	"github.com/pkg/errors"
)

var cmd = &cmdapp.Command{
	UsageLine: "p.len [-t|--tree <treefile>] <dataset>",
	Short:     "print the length of a tree",
	Long: `
Command p.len reads a tree in parenthetical format and prints its
length under parsimony.

The tree will be read from the standard input, unless the option
-t or --tree is defined with a tree file.

Options are:

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

func register(c *cmdapp.Command) {
	c.Flag.StringVar(&treefile, "tree", "", "")
	c.Flag.StringVar(&treefile, "t", "", "")
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

	m, err := matrix.NewMatrix(f)
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

	tr, err := parsimony.ReadTree(tf, m)
	if err != nil {
		return errors.Wrapf(err, "%s: when parsing tree", c.Name())
	}
	fmt.Printf("# Tree Length:\n%d\n", tr.Cost())
	return nil
}
