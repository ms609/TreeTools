# Identify consistent / conflicting splits

`SplitConsistent()` and `SplitConflict()` determine whether a series of
splits `haystack` are consistent with or contradict the focal split
`needle`.

## Usage

``` r
SplitConsistent(needle, haystack)

SplitConflicts(needle, haystack)
```

## Arguments

- needle:

  Splits object containing the single split to evaluate

- haystack:

  Splits object, or list thereof, containing the splits to compare
  against `needle`.

## Value

`SplitConsistent()` returns a list of logical vectors. Each list item
corresponds to an entry in `haystack`, reporting whether each split is
consistent with (`TRUE`) or in conflict with (`FALSE`) `needle`.
`SplitConflicts()` returns the inverse.

## See also

Other split manipulation functions:
[`DropTip()`](https://ms609.github.io/TreeTools/reference/DropTip.md),
[`Subsplit()`](https://ms609.github.io/TreeTools/reference/Subsplit.md),
[`TrivialSplits()`](https://ms609.github.io/TreeTools/reference/TrivialSplits.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
splits1 <- as.Splits(BalancedTree(8))
splits2 <- as.Splits(PectinateTree(8))
summary(splits1[[4]])
#> 1 bipartition split dividing 8 tips, t1 .. t8
#>      12345678
#>  14  ....**..
#> 
#>  Tip 1: t1    Tip 2: t2   Tip 3: t3   Tip 4: t4   Tip 5: t5  
#>  Tip 6: t6    Tip 7: t7   Tip 8: t8  
SplitConsistent(splits1[[4]], splits2)
#> [[1]]
#> [1]  TRUE  TRUE  TRUE FALSE  TRUE
#> 
SplitConflicts(splits1[[4]], list(splits1, splits2))
#> [[1]]
#> [1] FALSE FALSE FALSE FALSE FALSE
#> 
#> [[2]]
#> [1] FALSE FALSE FALSE  TRUE FALSE
#> 
```
