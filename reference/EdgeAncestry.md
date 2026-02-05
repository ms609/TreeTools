# Ancestors of an edge

Quickly identify edges that are "ancestral" to a particular edge in a
tree.

## Usage

``` r
EdgeAncestry(edge, parent, child, stopAt = (parent == min(parent)))
```

## Arguments

- edge:

  Integer specifying the number of the edge whose child edges should be
  returned.

- parent:

  Integer vector corresponding to the first column of the edge matrix of
  a tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html),
  i.e. `tree[["edge"]][, 1]`

- child:

  Integer vector corresponding to the second column of the edge matrix
  of a tree of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html), i.e.
  `tree[["edge"]][, 2]`.

- stopAt:

  Integer or logical vector specifying the edge(s) at which to terminate
  the search; defaults to the edges with the smallest parent, which will
  be the root edges if nodes are numbered
  [Cladewise](https://ms609.github.io/TreeTools/reference/Reorder.md) or
  in [Preorder](https://ms609.github.io/TreeTools/reference/Reorder.md).

## Value

`EdgeAncestry()` returns a logical vector stating whether each edge in
turn is a descendant of the specified edge.

## See also

Other tree navigation:
[`AncestorEdge()`](https://ms609.github.io/TreeTools/reference/AncestorEdge.md),
[`CladeSizes()`](https://ms609.github.io/TreeTools/reference/CladeSizes.md),
[`DescendantEdges()`](https://ms609.github.io/TreeTools/reference/DescendantEdges.md),
[`EdgeDistances()`](https://ms609.github.io/TreeTools/reference/EdgeDistances.md),
[`ListAncestors()`](https://ms609.github.io/TreeTools/reference/ListAncestors.md),
[`MRCA()`](https://ms609.github.io/TreeTools/reference/MRCA.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/reference/MatchEdges.md),
[`NDescendants()`](https://ms609.github.io/TreeTools/reference/NDescendants.md),
[`NodeDepth()`](https://ms609.github.io/TreeTools/reference/NodeDepth.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/reference/NodeNumbers.md),
[`NodeOrder()`](https://ms609.github.io/TreeTools/reference/NodeOrder.md),
[`RootNode()`](https://ms609.github.io/TreeTools/reference/RootNode.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tree <- PectinateTree(6)
plot(tree)
ape::edgelabels()

parent <- tree$edge[, 1]
child <- tree$edge[, 2]
EdgeAncestry(7, parent, child)
#>  [1] FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE
which(EdgeAncestry(7, parent, child, stopAt = 4))
#> [1] 4 6
```
