# Identify descendant edges

`DescendantEdges()` efficiently identifies edges that are "descended"
from edges in a tree.

`DescendantTips()` efficiently identifies leaves (external nodes) that
are "descended" from edges in a tree.

## Usage

``` r
DescendantEdges(
  parent,
  child,
  edge = NULL,
  node = NULL,
  nEdge = length(parent),
  includeSelf = TRUE
)

DescendantTips(parent, child, edge = NULL, node = NULL, nEdge = length(parent))
```

## Arguments

- parent:

  Integer vector corresponding to the first column of the edge matrix of
  a tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html),
  i.e. `tree[["edge"]][, 1]`

- child:

  Integer vector corresponding to the second column of the edge matrix
  of a tree of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html), i.e.
  `tree[["edge"]][, 2]`.

- edge:

  Integer specifying the number of the edge whose children are required
  (see [`edgelabels()`](https://rdrr.io/pkg/ape/man/nodelabels.html)).

- node:

  Integer specifying the number(s) of nodes whose children are required.
  Specify `0` to return all nodes. If `NULL` (the default), the `edge`
  parameter will be used instead.

- nEdge:

  number of edges (calculated from `length(parent)` if not supplied).

- includeSelf:

  Logical specifying whether to mark `edge` as its own descendant.

## Value

`DescendantEdges()` returns a logical vector stating whether each edge
in turn is the specified edge (if `includeSelf = TRUE`) or one of its
descendants.

`DescendantTips()` returns a logical vector stating whether each leaf in
turn is a descendant of the specified edge.

## See also

Other tree navigation:
[`AncestorEdge()`](https://ms609.github.io/TreeTools/dev/reference/AncestorEdge.md),
[`CladeSizes()`](https://ms609.github.io/TreeTools/dev/reference/CladeSizes.md),
[`EdgeAncestry()`](https://ms609.github.io/TreeTools/dev/reference/EdgeAncestry.md),
[`EdgeDistances()`](https://ms609.github.io/TreeTools/dev/reference/EdgeDistances.md),
[`ListAncestors()`](https://ms609.github.io/TreeTools/dev/reference/ListAncestors.md),
[`MRCA()`](https://ms609.github.io/TreeTools/dev/reference/MRCA.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/dev/reference/MatchEdges.md),
[`NDescendants()`](https://ms609.github.io/TreeTools/dev/reference/NDescendants.md),
[`NodeDepth()`](https://ms609.github.io/TreeTools/dev/reference/NodeDepth.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/dev/reference/NodeNumbers.md),
[`NodeOrder()`](https://ms609.github.io/TreeTools/dev/reference/NodeOrder.md),
[`RootNode()`](https://ms609.github.io/TreeTools/dev/reference/RootNode.md)

## Examples

``` r
tree <- as.phylo(0, 6)
plot(tree)
desc <- DescendantEdges(tree$edge[, 1], tree$edge[, 2], edge = 5)
which(desc)
#> [1] 5 6 7
ape::edgelabels(bg = 3 + desc)
tips <- DescendantTips(tree$edge[, 1], tree$edge[, 2], edge = 5)
which(tips)
#> [1] 2 6
tiplabels(bg = 3 + tips)
```
