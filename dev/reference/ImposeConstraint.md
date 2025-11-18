# Force a tree to match a constraint

Modify a tree such that it matches a specified constraint. This is at
present a somewhat crude implementation that attempts to retain much of
the structure of `tree` whilst guaranteeing compatibility with each
entry in `constraint`.

## Usage

``` r
ImposeConstraint(tree, constraint)

AddUnconstrained(constraint, toAdd, asPhyDat = TRUE)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- constraint:

  Either an object of class `phyDat`, in which case returned trees will
  be perfectly compatible with each character in `constraint`; or a tree
  of class `phylo`, in which each node in `constraint` will occur in the
  returned tree. See
  [vignette](https://ms609.github.io/TreeSearch/articles/tree-search.html)
  for further examples.

- toAdd:

  Character vector specifying taxa to add to constraint.

- asPhyDat:

  Logical: if `TRUE`, return a `phyDat` object; if `FALSE`, return a
  matrix.

## Value

`ImposeConstraint()` returns a tree of class `phylo`, consistent with
`constraint`.

## Functions

- `AddUnconstrained()`: Expand a constraint to include unconstrained
  taxa.

## See also

Other tree manipulation:
[`AddTip()`](https://ms609.github.io/TreeTools/dev/reference/AddTip.md),
[`CollapseNode()`](https://ms609.github.io/TreeTools/dev/reference/CollapseNode.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/dev/reference/ConsensusWithout.md),
[`DropTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md),
[`KeptPaths()`](https://ms609.github.io/TreeTools/dev/reference/KeptPaths.md),
[`KeptVerts()`](https://ms609.github.io/TreeTools/dev/reference/KeptVerts.md),
[`LeafLabelInterchange()`](https://ms609.github.io/TreeTools/dev/reference/LeafLabelInterchange.md),
[`MakeTreeBinary()`](https://ms609.github.io/TreeTools/dev/reference/MakeTreeBinary.md),
[`Renumber()`](https://ms609.github.io/TreeTools/dev/reference/Renumber.md),
[`RenumberTips()`](https://ms609.github.io/TreeTools/dev/reference/RenumberTips.md),
[`RenumberTree()`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md),
[`RootTree()`](https://ms609.github.io/TreeTools/dev/reference/RootTree.md),
[`SortTree()`](https://ms609.github.io/TreeTools/dev/reference/SortTree.md),
[`Subtree()`](https://ms609.github.io/TreeTools/dev/reference/Subtree.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/dev/reference/TipTimedTree.md),
[`TrivialTree`](https://ms609.github.io/TreeTools/dev/reference/TrivialTree.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tips <- letters[1:9]
tree <- as.phylo(1, 9, tips)
plot(tree)


constraint <- StringToPhyDat("0000?1111 000111111 0000??110", tips, FALSE)
plot(ImposeConstraint(tree, constraint))
```
