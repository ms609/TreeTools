# Count descendants for each node in a tree

`NDescendants()` counts the number of nodes (including leaves) directly
descended from each node in a tree.

## Usage

``` r
NDescendants(tree)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

## Value

`NDescendants()` returns an integer listing the number of direct
descendants (leaves or internal nodes) for each node in a tree.

## See also

Other tree navigation:
[`AncestorEdge()`](https://ms609.github.io/TreeTools/dev/reference/AncestorEdge.md),
[`CladeSizes()`](https://ms609.github.io/TreeTools/dev/reference/CladeSizes.md),
[`DescendantEdges()`](https://ms609.github.io/TreeTools/dev/reference/DescendantEdges.md),
[`EdgeAncestry()`](https://ms609.github.io/TreeTools/dev/reference/EdgeAncestry.md),
[`EdgeDistances()`](https://ms609.github.io/TreeTools/dev/reference/EdgeDistances.md),
[`ListAncestors()`](https://ms609.github.io/TreeTools/dev/reference/ListAncestors.md),
[`MRCA()`](https://ms609.github.io/TreeTools/dev/reference/MRCA.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/dev/reference/MatchEdges.md),
[`NodeDepth()`](https://ms609.github.io/TreeTools/dev/reference/NodeDepth.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/dev/reference/NodeNumbers.md),
[`NodeOrder()`](https://ms609.github.io/TreeTools/dev/reference/NodeOrder.md),
[`RootNode()`](https://ms609.github.io/TreeTools/dev/reference/RootNode.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tree <- CollapseNode(BalancedTree(8), 12:15)
NDescendants(tree)
#> [1] 5 3 2
plot(tree)
nodelabels(NDescendants(tree))

```
