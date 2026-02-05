# Distance of each node from tree exterior

`NodeDepth()` evaluates how "deep" each node is within a tree.

## Usage

``` r
NodeDepth(x, shortest = FALSE, includeTips = TRUE)
```

## Arguments

- x:

  A tree of class `phylo`, its `$edge` property, or a list thereof.

- shortest:

  Logical specifying whether to calculate the length of the shortest
  away-from-root path to a leaf. If `FALSE`, the length of the longest
  such route will be returned.

- includeTips:

  Logical specifying whether to include leaves (each of depth zero) in
  return value.

## Value

`NodeDepth()` returns an integer vector specifying the depth of each
external and internal node in `x`.

## Details

For a rooted tree, the depth of a node is the minimum (if
`shortest = TRUE`) or maximum (`shortest = FALSE`) number of edges that
must be traversed, moving away from the root, to reach a leaf.

Unrooted trees are treated as if a root node occurs in the "middle" of
the tree, meaning the position that will minimise the maximum node
depth.

## See also

[`ape::node.depth`](https://rdrr.io/pkg/ape/man/node.depth.html) returns
the number of tips descended from a node.

Other tree navigation:
[`AncestorEdge()`](https://ms609.github.io/TreeTools/reference/AncestorEdge.md),
[`CladeSizes()`](https://ms609.github.io/TreeTools/reference/CladeSizes.md),
[`DescendantEdges()`](https://ms609.github.io/TreeTools/reference/DescendantEdges.md),
[`EdgeAncestry()`](https://ms609.github.io/TreeTools/reference/EdgeAncestry.md),
[`EdgeDistances()`](https://ms609.github.io/TreeTools/reference/EdgeDistances.md),
[`ListAncestors()`](https://ms609.github.io/TreeTools/reference/ListAncestors.md),
[`MRCA()`](https://ms609.github.io/TreeTools/reference/MRCA.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/reference/MatchEdges.md),
[`NDescendants()`](https://ms609.github.io/TreeTools/reference/NDescendants.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/reference/NodeNumbers.md),
[`NodeOrder()`](https://ms609.github.io/TreeTools/reference/NodeOrder.md),
[`RootNode()`](https://ms609.github.io/TreeTools/reference/RootNode.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tree <- CollapseNode(BalancedTree(10), c(12:13, 19))
plot(tree)
nodelabels(NodeDepth(tree, includeTips = FALSE))


```
