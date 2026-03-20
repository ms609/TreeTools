# Number of trees

These functions return the number of rooted or unrooted binary trees
consistent with a given pattern of splits.

## Usage

``` r
NRooted(tips)

NUnrooted(tips)

NRooted64(tips)

NUnrooted64(tips)

LnUnrooted(tips)

LnUnrooted.int(tips)

Log2Unrooted(tips)

Log2Unrooted.int(tips)

LnRooted(tips)

LnRooted.int(tips)

Log2Rooted(tips)

Log2Rooted.int(tips)

LnUnrootedSplits(...)

Log2UnrootedSplits(...)

NUnrootedSplits(...)

LnUnrootedMult(...)

Log2UnrootedMult(...)

NUnrootedMult(...)
```

## Arguments

- tips:

  Integer specifying the number of leaves.

- ...:

  Integer vector, or series of integers, listing the number of leaves in
  each split.

## Details

Functions starting `N` return the number of rooted or unrooted trees.
Replace this initial `N` with `Ln` for the natural logarithm of this
number; or `Log2` for its base 2 logarithm.

Calculations follow Cavalli-Sforza and Edwards (1967) and Carter et al.
(1990) , Theorem 2.

## Functions

- `NUnrooted()`: Number of unrooted trees

- `NRooted64()`: Exact number of rooted trees as 64-bit integer (13 \<
  `nTip` \< 19)

- `NUnrooted64()`: Exact number of unrooted trees as 64-bit integer (14
  \< `nTip` \< 20)

- `LnUnrooted()`: Log Number of unrooted trees

- `LnUnrooted.int()`: Log Number of unrooted trees (as integer)

- `LnRooted()`: Log Number of rooted trees

- `LnRooted.int()`: Log Number of rooted trees (as integer)

- `NUnrootedSplits()`: Number of unrooted trees consistent with a
  bipartition split.

- `NUnrootedMult()`: Number of unrooted trees consistent with a
  multi-partition split.

## References

Carter M, Hendy M, Penny D, Székely LA, Wormald NC (1990). “On the
distribution of lengths of evolutionary trees.” *SIAM Journal on
Discrete Mathematics*, **3**(1), 38–47.
[doi:10.1137/0403005](https://doi.org/10.1137/0403005) .  
  
Cavalli-Sforza LL, Edwards AWF (1967). “Phylogenetic analysis: models
and estimation procedures.” *Evolution*, **21**(3), 550–570. ISSN
00143820,
[doi:10.1111/j.1558-5646.1967.tb03411.x](https://doi.org/10.1111/j.1558-5646.1967.tb03411.x)
.

## See also

Other tree information functions:
[`CladisticInfo()`](https://ms609.github.io/TreeTools/dev/reference/CladisticInfo.md),
[`TreesMatchingTree()`](https://ms609.github.io/TreeTools/dev/reference/TreesMatchingTree.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
NRooted(10)
#> [1] 34459425
NUnrooted(10)
#> [1] 2027025
LnRooted(10)
#> [1] 17.35529
LnUnrooted(10)
#> [1] 14.52208
Log2Unrooted(10)
#> [1] 20.95093
# Number of trees consistent with a character whose states are
# 00000 11111 222
NUnrootedMult(c(5,5,3))
#> [1] 694575

NUnrooted64(18)
#> integer64
#> [1] 191898783962510625
LnUnrootedSplits(c(2,4))
#> [1] 2.70805
LnUnrootedSplits(3, 3)
#> [1] 2.197225
Log2UnrootedSplits(c(2,4))
#> [1] 3.906891
Log2UnrootedSplits(3, 3)
#> [1] 3.169925
NUnrootedSplits(c(2,4))
#> [1] 15
NUnrootedSplits(3, 3)
#> [1] 9
```
