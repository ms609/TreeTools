# List ancestors

`ListAncestors()` reports all ancestors of a given node.

## Usage

``` r
ListAncestors(parent, child, node = NULL)

AllAncestors(parent, child)
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

- node:

  Integer giving the index of the node or tip whose ancestors are
  required, or `NULL` to return ancestors of all nodes.

## Value

If `node = NULL`, `ListAncestors()` returns a list. Each entry *i*
contains a vector containing, in order, the nodes encountered when
traversing the tree from node *i* to the root node. The last entry of
each member of the list is therefore the root node, with the exception
of the entry for the root node itself, which is a zero-length integer.

If `node` is an integer, `ListAncestors()` returns a vector of the
numbers of the nodes ancestral to the given `node`, including the root
node.

## Details

Note that if `node = NULL`, the tree's edges must be listed such that
each internal node (except the root) is listed as a child before it is
listed as a parent, i.e. its index in `child` is less than its index in
`parent`. This will be true of trees listed in
[Preorder](https://ms609.github.io/TreeTools/reference/Reorder.md).

## Functions

- `AllAncestors()`: Alias for `ListAncestors(node = NULL)`.

## See also

Implemented less efficiently in `phangorn:::Ancestors`, on which this
code is based.

Other tree navigation:
[`AncestorEdge()`](https://ms609.github.io/TreeTools/reference/AncestorEdge.md),
[`CladeSizes()`](https://ms609.github.io/TreeTools/reference/CladeSizes.md),
[`DescendantEdges()`](https://ms609.github.io/TreeTools/reference/DescendantEdges.md),
[`EdgeAncestry()`](https://ms609.github.io/TreeTools/reference/EdgeAncestry.md),
[`EdgeDistances()`](https://ms609.github.io/TreeTools/reference/EdgeDistances.md),
[`MRCA()`](https://ms609.github.io/TreeTools/reference/MRCA.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/reference/MatchEdges.md),
[`NDescendants()`](https://ms609.github.io/TreeTools/reference/NDescendants.md),
[`NodeDepth()`](https://ms609.github.io/TreeTools/reference/NodeDepth.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/reference/NodeNumbers.md),
[`NodeOrder()`](https://ms609.github.io/TreeTools/reference/NodeOrder.md),
[`RootNode()`](https://ms609.github.io/TreeTools/reference/RootNode.md)

Other tree navigation:
[`AncestorEdge()`](https://ms609.github.io/TreeTools/reference/AncestorEdge.md),
[`CladeSizes()`](https://ms609.github.io/TreeTools/reference/CladeSizes.md),
[`DescendantEdges()`](https://ms609.github.io/TreeTools/reference/DescendantEdges.md),
[`EdgeAncestry()`](https://ms609.github.io/TreeTools/reference/EdgeAncestry.md),
[`EdgeDistances()`](https://ms609.github.io/TreeTools/reference/EdgeDistances.md),
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
tree <- PectinateTree(5)
edge <- tree[["edge"]]

# Identify desired node with:
plot(tree)
nodelabels()
tiplabels()


# Ancestors of specific nodes:
ListAncestors(edge[, 1], edge[, 2], 4L)
#> [1] 9 8 7 6
ListAncestors(edge[, 1], edge[, 2], 8L)
#> [1] 7 6

# Ancestors of each node, if tree numbering system is uncertain:
lapply(seq_len(max(edge)), ListAncestors,
       parent = edge[, 1], child = edge[, 2])
#> [[1]]
#> [1] 6
#> 
#> [[2]]
#> [1] 7 6
#> 
#> [[3]]
#> [1] 8 7 6
#> 
#> [[4]]
#> [1] 9 8 7 6
#> 
#> [[5]]
#> [1] 9 8 7 6
#> 
#> [[6]]
#> integer(0)
#> 
#> [[7]]
#> [1] 6
#> 
#> [[8]]
#> [1] 7 6
#> 
#> [[9]]
#> [1] 8 7 6
#> 

# Ancestors of each node, if tree is in preorder:
ListAncestors(edge[, 1], edge[, 2])
#> [[1]]
#> [1] 6
#> 
#> [[2]]
#> [1] 7 6
#> 
#> [[3]]
#> [1] 8 7 6
#> 
#> [[4]]
#> [1] 9 8 7 6
#> 
#> [[5]]
#> [1] 9 8 7 6
#> 
#> [[6]]
#> integer(0)
#> 
#> [[7]]
#> [1] 6
#> 
#> [[8]]
#> [1] 7 6
#> 
#> [[9]]
#> [1] 8 7 6
#> 

# Alias:
AllAncestors(edge[, 1], edge[, 2])
#> [[1]]
#> [1] 6
#> 
#> [[2]]
#> [1] 7 6
#> 
#> [[3]]
#> [1] 8 7 6
#> 
#> [[4]]
#> [1] 9 8 7 6
#> 
#> [[5]]
#> [1] 9 8 7 6
#> 
#> [[6]]
#> integer(0)
#> 
#> [[7]]
#> [1] 6
#> 
#> [[8]]
#> [1] 7 6
#> 
#> [[9]]
#> [1] 8 7 6
#> 
```
