# Distance between edges

Number of nodes that must be traversed to navigate from each edge to
each other edge within a tree

## Usage

``` r
EdgeDistances(tree)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

## Value

`EdgeDistances()` returns a symmetrical matrix listing the number of
edges that must be traversed to travel from each numbered edge to each
other. The two edges straddling the root of a rooted tree are treated as
a single edge. Add a "root" tip using
[`AddTip()`](https://ms609.github.io/TreeTools/dev/reference/AddTip.md)
if the position of the root is significant.

## See also

Other tree navigation:
[`AncestorEdge()`](https://ms609.github.io/TreeTools/dev/reference/AncestorEdge.md),
[`CladeSizes()`](https://ms609.github.io/TreeTools/dev/reference/CladeSizes.md),
[`DescendantEdges()`](https://ms609.github.io/TreeTools/dev/reference/DescendantEdges.md),
[`EdgeAncestry()`](https://ms609.github.io/TreeTools/dev/reference/EdgeAncestry.md),
[`ListAncestors()`](https://ms609.github.io/TreeTools/dev/reference/ListAncestors.md),
[`MRCA()`](https://ms609.github.io/TreeTools/dev/reference/MRCA.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/dev/reference/MatchEdges.md),
[`NDescendants()`](https://ms609.github.io/TreeTools/dev/reference/NDescendants.md),
[`NodeDepth()`](https://ms609.github.io/TreeTools/dev/reference/NodeDepth.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/dev/reference/NodeNumbers.md),
[`NodeOrder()`](https://ms609.github.io/TreeTools/dev/reference/NodeOrder.md),
[`RootNode()`](https://ms609.github.io/TreeTools/dev/reference/RootNode.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tree <- BalancedTree(5)
plot(tree)
ape::edgelabels()


EdgeDistances(tree)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,]    0    1    2    2    1    0    1    1
#> [2,]    1    0    1    1    1    1    2    2
#> [3,]    2    1    0    1    2    2    3    3
#> [4,]    2    1    1    0    2    2    3    3
#> [5,]    1    1    2    2    0    1    2    2
#> [6,]    0    1    2    2    1    0    1    1
#> [7,]    1    2    3    3    2    1    0    1
#> [8,]    1    2    3    3    2    1    1    0
```
