# Paths present in reduced tree

Lists which paths present in a master tree are present when leaves are
dropped.

## Usage

``` r
KeptPaths(paths, keptVerts, all = TRUE)

# S3 method for class 'data.frame'
KeptPaths(paths, keptVerts, all = TRUE)

# S3 method for class 'matrix'
KeptPaths(paths, keptVerts, all = TRUE)
```

## Arguments

- paths:

  `data.frame` of paths in master tree, perhaps generated using
  [`PathLengths()`](https://ms609.github.io/TreeTools/reference/PathLengths.md).

- keptVerts:

  Logical specifying whether each entry is retained in the reduced tree,
  perhaps generated using
  [`KeptVerts()`](https://ms609.github.io/TreeTools/reference/KeptVerts.md).

- all:

  Logical: if `TRUE`, return all paths that occur in the reduced tree;
  if `FALSE`, return only those paths that correspond to a single edge.
  that correspond to edges in the reduced tree. Ignored if `paths` is a
  matrix.

## Value

`KeptPaths()` returns a logical vector specifying whether each path in
`paths` occurs when `keptVerts` vertices are retained.

## See also

Other tree manipulation:
[`AddTip()`](https://ms609.github.io/TreeTools/reference/AddTip.md),
[`CollapseNode()`](https://ms609.github.io/TreeTools/reference/CollapseNode.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/reference/ConsensusWithout.md),
[`DropTip()`](https://ms609.github.io/TreeTools/reference/DropTip.md),
[`ImposeConstraint()`](https://ms609.github.io/TreeTools/reference/ImposeConstraint.md),
[`KeptVerts()`](https://ms609.github.io/TreeTools/reference/KeptVerts.md),
[`LeafLabelInterchange()`](https://ms609.github.io/TreeTools/reference/LeafLabelInterchange.md),
[`MakeTreeBinary()`](https://ms609.github.io/TreeTools/reference/MakeTreeBinary.md),
[`Renumber()`](https://ms609.github.io/TreeTools/reference/Renumber.md),
[`RenumberTips()`](https://ms609.github.io/TreeTools/reference/RenumberTips.md),
[`RenumberTree()`](https://ms609.github.io/TreeTools/reference/Reorder.md),
[`RootTree()`](https://ms609.github.io/TreeTools/reference/RootTree.md),
[`SortTree()`](https://ms609.github.io/TreeTools/reference/SortTree.md),
[`Subtree()`](https://ms609.github.io/TreeTools/reference/Subtree.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/reference/TipTimedTree.md),
[`TrivialTree`](https://ms609.github.io/TreeTools/reference/TrivialTree.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
master <- BalancedTree(9)
paths <- PathLengths(master)
keptTips <- c(1, 5, 7, 9)
keptVerts <- KeptVerts(master, keptTips)
KeptPaths(paths, keptVerts)
#>  [1]  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [13] FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE
#> [25] FALSE FALSE  TRUE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
#> [37] FALSE  TRUE FALSE FALSE FALSE FALSE
paths[KeptPaths(paths, keptVerts, all = FALSE), ]
#>    start end length
#> 2     11   1      3
#> 16    11   5      2
#> 22    15   7      2
#> 28    15   9      2
#> 30    10  11      1
#> 38    10  15      1
```
