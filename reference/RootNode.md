# Which node is a tree's root?

`RootNode()` identifies the root node of a (rooted or unrooted)
phylogenetic tree. Unrooted trees are represented internally by a rooted
tree with a polytomy at the root.

## Usage

``` r
RootNode(x)
```

## Arguments

- x:

  A tree of class `phylo`, or its edge matrix; or a list or `multiPhylo`
  object containing multiple trees.

## Value

`RootNode()` returns an integer denoting the root node for each tree.
Badly conformed trees trigger an error.

## See also

Test whether a tree is rooted:
[`TreeIsRooted()`](https://ms609.github.io/TreeTools/reference/TreeIsRooted.md)

[`phangorn::getRoot()`](https://klausvigo.github.io/phangorn/reference/midpoint.html)

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
[`NodeDepth()`](https://ms609.github.io/TreeTools/reference/NodeDepth.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/reference/NodeNumbers.md),
[`NodeOrder()`](https://ms609.github.io/TreeTools/reference/NodeOrder.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
RootNode(BalancedTree(8))
#> [1] 9
RootNode(UnrootTree(BalancedTree(8)))
#> [1] 9

```
