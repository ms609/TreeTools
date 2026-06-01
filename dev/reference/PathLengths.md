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
#> 1      7   1 2.1695031
#> 2      8   1 1.5267076
#> 3      9   1 0.5980924
#> 4      7   2 2.6583392
#> 5      8   2 2.0155437
#> 6      9   2 1.0869285
#> 7     10   2 0.5260277
#> 8      7   3 3.1174067
#> 9      8   3 2.4746112
#> 10     9   3 1.5459960
#> 11    10   3 0.9850952
#> 12     7   4 1.1504373
#> 13     8   4 0.5076418
#> 14     7   5 1.2843293
#> 15    11   5 0.6015412
#> 16     7   6 0.9216568
#> 17    11   6 0.2388687
#> 18     7   8 0.6427955
#> 19     7   9 1.5714107
#> 20     8   9 0.9286152
#> 21     7  10 2.1323114
#> 22     8  10 1.4895159
#> 23     9  10 0.5609007
#> 24     7  11 0.6827881
```
