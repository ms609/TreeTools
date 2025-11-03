# Number of trees containing a tree

`TreesMatchingTree()` calculates the number of unrooted binary trees
that are consistent with a tree topology on the same leaves.

## Usage

``` r
TreesMatchingTree(tree)

LnTreesMatchingTree(tree)

Log2TreesMatchingTree(tree)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

## Value

`TreesMatchingTree()` returns a numeric specifying the number of
unrooted binary trees that contain all the edges present in the input
tree.

`LnTreesMatchingTree()` gives the natural logarithm of this number.

## Details

Remember to unroot a tree first if the position of its root is
arbitrary.

## See also

Other tree information functions:
[`CladisticInfo()`](https://ms609.github.io/TreeTools/reference/CladisticInfo.md),
[`NRooted()`](https://ms609.github.io/TreeTools/reference/NRooted.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
partiallyResolvedTree <- CollapseNode(BalancedTree(8), 12:15)
TreesMatchingTree(partiallyResolvedTree)
#> [1] 45
LnTreesMatchingTree(partiallyResolvedTree)
#> [1] 3.806662

# Number of rooted trees:
rootedTree <- AddTip(partiallyResolvedTree, where = 0)
TreesMatchingTree(partiallyResolvedTree)
#> [1] 45
```
