# Reorder edges of a phylogenetic tree

Wrappers for the C functions called by
`ape::`[`reorder.phylo`](https://rdrr.io/pkg/ape/man/reorder.phylo.html).
These call the C functions directly, so are faster – but don't perform
as many checks on user input. Bad input could crash R.

## Usage

``` r
NeworderPruningwise(nTip, nNode, parent, child, nEdge)

NeworderPhylo(nTip, parent, child, nEdge, whichwise)
```

## Arguments

- nTip, nNode, nEdge:

  Integer specifying the number of tips, nodes and edges in the input
  tree.

- parent:

  Integer vector corresponding to the first column of the edge matrix of
  a tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html),
  i.e. `tree[["edge"]][, 1]`

- child:

  Integer vector corresponding to the second column of the edge matrix
  of a tree of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html), i.e.
  `tree[["edge"]][, 2]`.

- whichwise:

  Integer specifying whether to order edges (1) cladewise; or (2) in
  postorder.

## Value

`NeworderPruningwise` returns an integer vector specifying the
pruningwise order of edges within a tree.

`NeworderPhylo` returns an integer vector specifying the order of edges
under the ordering sequence specified by `whichwise`.

## See also

Other C wrappers:
[`RenumberTree()`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md)

## Author

- C algorithm: Emmanuel Paradis

- R wrapper: Martin R. Smith

## Examples

``` r
nTip <- 8L
tree <- BalancedTree(nTip)
edge <- tree[["edge"]]
pruningwise <- NeworderPruningwise(nTip, tree$Nnode, edge[, 1], edge[, 2],
                                   dim(edge)[1])
cladewise <- NeworderPhylo(nTip, edge[, 1], edge[, 2], dim(edge)[1], 1L)
postorder <- NeworderPhylo(nTip, edge[, 1], edge[, 2], dim(edge)[1], 2L)

tree[["edge"]] <- tree[["edge"]][pruningwise, ]
```
