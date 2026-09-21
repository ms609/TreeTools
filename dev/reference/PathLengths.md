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
  [`Preorder`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md).

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
[`Cherries()`](https://ms609.github.io/TreeTools/dev/reference/Cherries.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/dev/reference/ConsensusWithout.md),
[`EdgeRatio()`](https://ms609.github.io/TreeTools/dev/reference/EdgeRatio.md),
[`LongBranch()`](https://ms609.github.io/TreeTools/dev/reference/LongBranch.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/dev/reference/MatchEdges.md),
[`NSplits()`](https://ms609.github.io/TreeTools/dev/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/dev/reference/NTip.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/dev/reference/NodeNumbers.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/dev/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/dev/reference/TipLabels.md),
[`TreeIsRooted()`](https://ms609.github.io/TreeTools/dev/reference/TreeIsRooted.md),
[`Treeness()`](https://ms609.github.io/TreeTools/dev/reference/Treeness.md)

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
#>    start end    length
#> 1      7   1 1.3248956
#> 2      8   1 0.7617968
#> 3      7   2 2.0162235
#> 4      8   2 1.4531247
#> 5      9   2 0.7078273
#> 6      7   3 2.2952045
#> 7      8   3 1.7321056
#> 8      9   3 0.9868082
#> 9      7   4 0.8084836
#> 10    10   4 0.3568527
#> 11    11   4 0.2386353
#> 12     7   5 1.0403711
#> 13    10   5 0.5887402
#> 14    11   5 0.4705228
#> 15     7   6 0.8039540
#> 16    10   6 0.3523231
#> 17     7   8 0.5630988
#> 18     7   9 1.3083962
#> 19     8   9 0.7452974
#> 20     7  10 0.4516309
#> 21     7  11 0.5698483
#> 22    10  11 0.1182174
```
