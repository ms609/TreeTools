# Number of edges incident to each node in a tree

`NodeOrder()` calculates the order of each node: the number of edges
incident to it in a tree. This value includes the root edge in rooted
trees.

## Usage

``` r
NodeOrder(x, includeAncestor = TRUE, internalOnly = FALSE)
```

## Arguments

- x:

  A tree of class `phylo`, its `$edge` property, or a list thereof.

- includeAncestor:

  Logical specifying whether to count edge leading to ancestral node in
  calculation of order.

- internalOnly:

  Logical specifying whether to restrict to results to internal nodes,
  i.e. to omit leaves. Irrelevant if `includeAncestor = FALSE`.

## Value

`NodeOrder()` returns an integer listing the order of each node; entries
are named with the number of each node.

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
[`NDescendants()`](https://ms609.github.io/TreeTools/dev/reference/NDescendants.md),
[`NodeDepth()`](https://ms609.github.io/TreeTools/dev/reference/NodeDepth.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/dev/reference/NodeNumbers.md),
[`RootNode()`](https://ms609.github.io/TreeTools/dev/reference/RootNode.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tree <- CollapseNode(BalancedTree(8), 12:15)
NodeOrder(tree)
#>  [1] 1 1 1 1 1 1 1 1 5 4 3
plot(tree)
nodelabels(NodeOrder(tree, internalOnly = TRUE))

```
