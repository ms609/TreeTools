# Most recent common ancestor

`MRCA()` calculates the last common ancestor of specified nodes.

## Usage

``` r
MRCA(x1, x2, ancestors)
```

## Arguments

- x1, x2:

  Integer specifying index of leaves or nodes whose most recent common
  ancestor should be found.

- ancestors:

  List of ancestors for each node in a tree. Perhaps produced by
  [`ListAncestors()`](https://ms609.github.io/TreeTools/reference/ListAncestors.md).

## Value

`MRCA()` returns an integer specifying the node number of the last
common ancestor of `x1` and `x2`.

## Details

`MRCA()` requires that node values within a tree increase away from the
root, which will be true of trees listed in `Preorder`. No warnings will
be given if trees do not fulfil this requirement.

## See also

Other tree navigation:
[`AncestorEdge()`](https://ms609.github.io/TreeTools/reference/AncestorEdge.md),
[`CladeSizes()`](https://ms609.github.io/TreeTools/reference/CladeSizes.md),
[`DescendantEdges()`](https://ms609.github.io/TreeTools/reference/DescendantEdges.md),
[`EdgeAncestry()`](https://ms609.github.io/TreeTools/reference/EdgeAncestry.md),
[`EdgeDistances()`](https://ms609.github.io/TreeTools/reference/EdgeDistances.md),
[`ListAncestors()`](https://ms609.github.io/TreeTools/reference/ListAncestors.md),
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
tree <- BalancedTree(7)

# Verify that node numbering increases away from root
plot(tree)
nodelabels()


# ListAncestors expects a tree in Preorder
tree <- Preorder(tree)
edge <- tree$edge
ancestors <- ListAncestors(edge[, 1], edge[, 2])
MRCA(1, 4, ancestors)
#> [1] 9

# If a tree must be in postorder, use:
tree <- Postorder(tree)
edge <- tree$edge
ancestors <- lapply(seq_len(max(edge)), ListAncestors,
                    parent = edge[, 1], child = edge[, 2])
```
