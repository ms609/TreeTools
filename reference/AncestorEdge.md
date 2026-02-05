# Ancestral edge

Ancestral edge

## Usage

``` r
AncestorEdge(edge, parent, child)
```

## Arguments

- edge:

  Number of an edge

- parent:

  Integer vector corresponding to the first column of the edge matrix of
  a tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html),
  i.e. `tree[["edge"]][, 1]`

- child:

  Integer vector corresponding to the second column of the edge matrix
  of a tree of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html), i.e.
  `tree[["edge"]][, 2]`.

## Value

`AncestorEdge` returns a logical vector identifying whether each edge is
the immediate ancestor of the given edge.

## See also

Other tree navigation:
[`CladeSizes()`](https://ms609.github.io/TreeTools/reference/CladeSizes.md),
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
parent <- tree$edge[, 1]
child <- tree$edge[, 2]
plot(tree)
ape::edgelabels()

AncestorEdge(5, parent, child)
#>  [1]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
which(AncestorEdge(5, parent, child))
#> [1] 1
```
