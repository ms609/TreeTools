# Match nodes and edges between trees

`MatchNodes()` and `MatchEdges()` matches nodes or edges in one tree to
entries in the second that denote a clade with identical tip labels.

## Usage

``` r
MatchEdges(x, table, nomatch = NA_integer_)

MatchNodes(x, table, nomatch = NA_integer_, tips = FALSE)
```

## Arguments

- x:

  Tree whose nodes are to be matched.

- table:

  Tree containing nodes to be matched against.

- nomatch:

  Integer value that will be used in place of `NA` in the case where no
  match is found.

- tips:

  Logical specifying whether to return matches for tips; unless `TRUE`,
  only the matches for internal nodes will be returned.

## Value

`MatchNodes()` returns an integer vector giving, for each internal node
of `x` (and, if `tips = TRUE`, each tip, listed first), the number of
the node in `table` that subtends the same clade, or `nomatch` if no
node corresponds.

`MatchEdges()` returns an integer vector giving, for each edge of `x`,
the index of the edge in `table` that subtends the same clade, or
`nomatch`. Where several nodes or edges of `table` denote the same clade
(because of nodes with a single child), the first is returned.

## Details

`MatchNodes()` treats clades as identical only if they contain exactly
the same tips. `MatchEdges()` disregards tips of `table` that are absent
from `x`.

## References

There are no references for Rd macro `\insertAllCites` on this help
page.

## See also

Other tree navigation:
[`AncestorEdge()`](https://ms609.github.io/TreeTools/reference/AncestorEdge.md),
[`CladeSizes()`](https://ms609.github.io/TreeTools/reference/CladeSizes.md),
[`DescendantEdges()`](https://ms609.github.io/TreeTools/reference/DescendantEdges.md),
[`EdgeAncestry()`](https://ms609.github.io/TreeTools/reference/EdgeAncestry.md),
[`EdgeDistances()`](https://ms609.github.io/TreeTools/reference/EdgeDistances.md),
[`ListAncestors()`](https://ms609.github.io/TreeTools/reference/ListAncestors.md),
[`MRCA()`](https://ms609.github.io/TreeTools/reference/MRCA.md),
[`NDescendants()`](https://ms609.github.io/TreeTools/reference/NDescendants.md),
[`NodeDepth()`](https://ms609.github.io/TreeTools/reference/NodeDepth.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/reference/NodeNumbers.md),
[`NodeOrder()`](https://ms609.github.io/TreeTools/reference/NodeOrder.md),
[`PaintTree()`](https://ms609.github.io/TreeTools/reference/PaintTree.md),
[`RootNode()`](https://ms609.github.io/TreeTools/reference/RootNode.md)

Other tree properties:
[`Cherries()`](https://ms609.github.io/TreeTools/reference/Cherries.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/reference/ConsensusWithout.md),
[`EdgeRatio()`](https://ms609.github.io/TreeTools/reference/EdgeRatio.md),
[`LongBranch()`](https://ms609.github.io/TreeTools/reference/LongBranch.md),
[`NSplits()`](https://ms609.github.io/TreeTools/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/reference/NTip.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/reference/NodeNumbers.md),
[`PathLengths()`](https://ms609.github.io/TreeTools/reference/PathLengths.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/reference/TipLabels.md),
[`TreeIsRooted()`](https://ms609.github.io/TreeTools/reference/TreeIsRooted.md),
[`Treeness()`](https://ms609.github.io/TreeTools/reference/Treeness.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
MatchNodes(BalancedTree(8), RootTree(BalancedTree(8)))
#> [1]  9 10 11 12 13 14 15
MatchEdges(BalancedTree(8), RootTree(BalancedTree(8)))
#>  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14
```
