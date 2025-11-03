# Clade sizes

`CladeSizes()` reports the number of nodes in each clade in a tree.

## Usage

``` r
CladeSizes(tree, internal = FALSE, nodes = NULL)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- internal:

  Logical specifying whether internal nodes should be counted towards
  the size of each clade.

- nodes:

  Integer specifying indices of nodes at the base of clades whose sizes
  should be returned. If unspecified, counts will be provided for all
  nodes (including leaves).

## Value

`CladeSizes()` returns the number of nodes (including leaves) that are
descended from each node, not including the node itself.

## See also

Other tree navigation:
[`AncestorEdge()`](https://ms609.github.io/TreeTools/reference/AncestorEdge.md),
[`DescendantEdges()`](https://ms609.github.io/TreeTools/reference/DescendantEdges.md),
[`EdgeAncestry()`](https://ms609.github.io/TreeTools/reference/EdgeAncestry.md),
[`EdgeDistances()`](https://ms609.github.io/TreeTools/reference/EdgeDistances.md),
[`ListAncestors()`](https://ms609.github.io/TreeTools/reference/ListAncestors.md),
[`MRCA()`](https://ms609.github.io/TreeTools/reference/MRCA.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/reference/MatchEdges.md),
[`NDescendants()`](https://ms609.github.io/TreeTools/reference/NDescendants.md),
[`NodeDepth()`](https://ms609.github.io/TreeTools/reference/NodeDepth.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/reference/NodeNumbers.md),
[`NodeOrder()`](https://ms609.github.io/TreeTools/reference/NodeOrder.md),
[`RootNode()`](https://ms609.github.io/TreeTools/reference/RootNode.md)

## Examples

``` r
tree <- BalancedTree(6)
plot(tree)
ape::nodelabels()

CladeSizes(tree, nodes = c(1, 8, 9))
#> [1] 1 3 2
```
