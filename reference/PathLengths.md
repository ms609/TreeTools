# Calculate length of paths between each pair of vertices within tree

Given a weighted rooted tree `tree`, `PathLengths()` returns the
distance from each vertex to each of its descendant vertices.

## Usage

``` r
PathLengths(tree, fullMatrix = FALSE, use.na = TRUE)
```

## Arguments

- tree:

  Original tree of class `phylo`, in
  [`Preorder`](https://ms609.github.io/TreeTools/reference/Reorder.md).

- fullMatrix:

  Logical specifying return format; see "value" section\`.

- use.na:

  Logical specifying whether to set non-existent paths to `NA`, or to
  leave uninitialized. Set to `FALSE` to maximize performance.

## Value

If `fullMatrix = TRUE`, `PathLengths()` returns a square matrix in which
entry `[i, j]` denotes the distance from internal node `i` to the
descendant vertex `j`. Vertex pairs without a continuous directed path
are denoted `NA` if `use.na` is `TRUE`. If `fullMatrix = FALSE`,
`PathLengths()` returns a `data.frame` with three columns: `start` lists
the deepest node in each path (i.e. that closest to the root); `end`
lists the shallowest node (i.e. that closest to a leaf); `length` lists
the total length of that path.

## See also

Other tree properties:
[`Cherries()`](https://ms609.github.io/TreeTools/reference/Cherries.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/reference/ConsensusWithout.md),
[`EdgeRatio()`](https://ms609.github.io/TreeTools/reference/EdgeRatio.md),
[`LongBranch()`](https://ms609.github.io/TreeTools/reference/LongBranch.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/reference/MatchEdges.md),
[`NSplits()`](https://ms609.github.io/TreeTools/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/reference/NTip.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/reference/NodeNumbers.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/reference/TipLabels.md),
[`TreeIsRooted()`](https://ms609.github.io/TreeTools/reference/TreeIsRooted.md),
[`Treeness()`](https://ms609.github.io/TreeTools/reference/Treeness.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tree <- rtree(6)
plot(tree)
add.scale.bar()
nodelabels()
tiplabels()

PathLengths(tree)
#>    start end     length
#> 1      7   1 1.41261280
#> 2      8   1 0.92826328
#> 3      9   1 0.75482094
#> 4      7   2 1.11168735
#> 5      8   2 0.62733782
#> 6      9   2 0.45389549
#> 7      7   3 0.99551931
#> 8      8   3 0.51116978
#> 9      7   4 1.03191525
#> 10    10   4 0.82437014
#> 11    11   4 0.59571200
#> 12     7   5 1.01107545
#> 13    10   5 0.80353034
#> 14    11   5 0.57487220
#> 15     7   6 0.28460949
#> 16    10   6 0.07706438
#> 17     7   8 0.48434952
#> 18     7   9 0.65779186
#> 19     8   9 0.17344233
#> 20     7  10 0.20754511
#> 21     7  11 0.43620326
#> 22    10  11 0.22865814
```
