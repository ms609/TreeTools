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

## Details

The current implementation is potentially inefficient. Please contact
the maintainer to request a more efficient implementation if this
function is proving a bottleneck.

## See also

Other tree navigation:
[`AncestorEdge()`](https://ms609.github.io/TreeTools/dev/reference/AncestorEdge.md),
[`CladeSizes()`](https://ms609.github.io/TreeTools/dev/reference/CladeSizes.md),
[`DescendantEdges()`](https://ms609.github.io/TreeTools/dev/reference/DescendantEdges.md),
[`EdgeAncestry()`](https://ms609.github.io/TreeTools/dev/reference/EdgeAncestry.md),
[`EdgeDistances()`](https://ms609.github.io/TreeTools/dev/reference/EdgeDistances.md),
[`ListAncestors()`](https://ms609.github.io/TreeTools/dev/reference/ListAncestors.md),
[`MRCA()`](https://ms609.github.io/TreeTools/dev/reference/MRCA.md),
[`NDescendants()`](https://ms609.github.io/TreeTools/dev/reference/NDescendants.md),
[`NodeDepth()`](https://ms609.github.io/TreeTools/dev/reference/NodeDepth.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/dev/reference/NodeNumbers.md),
[`NodeOrder()`](https://ms609.github.io/TreeTools/dev/reference/NodeOrder.md),
[`RootNode()`](https://ms609.github.io/TreeTools/dev/reference/RootNode.md)

Other tree properties:
[`Cherries()`](https://ms609.github.io/TreeTools/dev/reference/Cherries.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/dev/reference/ConsensusWithout.md),
[`LongBranch()`](https://ms609.github.io/TreeTools/dev/reference/LongBranch.md),
[`NSplits()`](https://ms609.github.io/TreeTools/dev/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/dev/reference/NTip.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/dev/reference/NodeNumbers.md),
[`PathLengths()`](https://ms609.github.io/TreeTools/dev/reference/PathLengths.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/dev/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/dev/reference/TipLabels.md),
[`TreeIsRooted()`](https://ms609.github.io/TreeTools/dev/reference/TreeIsRooted.md),
[`Treeness()`](https://ms609.github.io/TreeTools/dev/reference/Treeness.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
MatchNodes(BalancedTree(8), RootTree(BalancedTree(8)))
#> [1]  9 10 11 12 13 14 15
```
