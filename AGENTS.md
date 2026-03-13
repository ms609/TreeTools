TreeTools underpins a hierarchy of packages, all with the same maintainer;
these draw in C++ headers and many exported R functions.

These are likely checked out locally in sister directories to this one.

## Main stack:
PlotTools
TreeTools (foundation-level)
TreeDist
TreeSearch ("top of the stack")

## Auxiliary packages:
Quartet
Rogue

## Complementary packages:
Ternary

Each package contains CONTRIBUTING.md files that detail code style conventions.

## Common data structures

Trees are represented ape's as.phylo, with edges listed as a two-column matrix
(parent node ID, child ID).
My Preorder ordering guarantees a particular sequence of edges and numbering
of internal nodes for any topologically identical tree.

Splits objects are defined in `as.Splits()`, and denote split membership as
binary 0/1 in an underlying `raw` object.

## Workflow requirements

- After completing each optimization or user-visible change, update `NEWS.md`
  before moving on to the next task.
- Increment the `.900X` dev version suffix in `DESCRIPTION` with each
  `NEWS.md` update.
- Check that existing tests cover all new code. (The GHA test suite uses codecov.)
