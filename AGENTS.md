# NA

> **Before starting work in this directory, read
> [`../AGENTS.md`](https://ms609.github.io/AGENTS.md)** for multi-agent
> coordination rules, build/test infrastructure, GHA workflows, and
> worktree discipline. That file is the authoritative reference for all
> cross-package agent operations.

TreeTools underpins a hierarchy of packages, all with the same
maintainer; these draw in C++ headers and many exported R functions.

These are likely checked out locally in sister directories to this one.

## Main stack:

PlotTools TreeTools (foundation-level) TreeDist TreeSearch (“top of the
stack”)

## Auxiliary packages:

Quartet Rogue

## Complementary packages:

Ternary

Each package contains CONTRIBUTING.md files that detail code style
conventions.

## Common data structures

Trees are represented ape’s as.phylo, with edges listed as a two-column
matrix (parent node ID, child ID). My Preorder ordering guarantees a
particular sequence of edges and numbering of internal nodes for any
topologically identical tree.

Splits objects are defined in
[`as.Splits()`](https://ms609.github.io/TreeTools/reference/Splits.md),
and denote split membership as binary 0/1 in an underlying `raw` object.

## Workflow requirements

- After completing each optimization or user-visible change, update
  `NEWS.md` before moving on to the next task.
- Increment the `.900X` dev version suffix in `DESCRIPTION` with each
  `NEWS.md` update.
- All new and changed code must have test coverage. The GHA test suite
  uses codecov; uncovered lines will block the PR. Cover happy paths,
  error branches, and edge cases (e.g. early returns).

## Optimization notes

- `descendant_edges_single()`: the O(n_edge) linear scan per node looks
  theoretically O(n²), but benchmarking (CSR index, vector-of-vectors)
  showed the original is faster at all practical sizes (up to 50k tips)
  due to cache-friendly sequential access over contiguous Rcpp memory.
  Not worth optimizing.
