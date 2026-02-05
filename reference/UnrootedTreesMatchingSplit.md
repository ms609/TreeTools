# Number of trees consistent with split

Calculates the number of unrooted bifurcating trees consistent with the
specified multi-partition split, using theorem two of Carter et al.
(1990) .

## Usage

``` r
UnrootedTreesMatchingSplit(...)

LnUnrootedTreesMatchingSplit(...)

Log2UnrootedTreesMatchingSplit(...)
```

## Arguments

- ...:

  A series or vector of integers listing the number of tips in each of a
  number of tree splits (e.g. bipartitions). For example, `3, 5` states
  that a character divides a set of eight tips into a group of three and
  a group of five.

## Value

`UnrootedTreesMatchingSplit()` returns an integer specifying the number
of unrooted bifurcating trees consistent with the specified split.

## References

Carter M, Hendy M, Penny D, Székely LA, Wormald NC (1990). “On the
distribution of lengths of evolutionary trees.” *SIAM Journal on
Discrete Mathematics*, **3**(1), 38–47.
[doi:10.1137/0403005](https://doi.org/10.1137/0403005) .

## See also

Other split information functions:
[`CharacterInformation()`](https://ms609.github.io/TreeTools/reference/CharacterInformation.md),
[`SplitInformation()`](https://ms609.github.io/TreeTools/reference/SplitInformation.md),
[`SplitMatchProbability()`](https://ms609.github.io/TreeTools/reference/SplitMatchProbability.md),
[`TreesMatchingSplit()`](https://ms609.github.io/TreeTools/reference/TreesMatchingSplit.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
UnrootedTreesMatchingSplit(c(3, 5))
#> [1] 315
UnrootedTreesMatchingSplit(3, 2, 1, 2)
#> [1] 297
```
