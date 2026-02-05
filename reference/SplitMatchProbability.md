# Probability of matching this well

(`Ln`)`SplitMatchProbability()`calculates the probability that two
random splits of the sizes provided will be at least as similar as the
two specified.

## Usage

``` r
SplitMatchProbability(split1, split2)

LnSplitMatchProbability(split1, split2)
```

## Arguments

- split1, split2:

  Logical vectors listing terminals in same order, such that each
  terminal is identified as a member of the ingroup (`TRUE`) or outgroup
  (`FALSE`) of the respective bipartition split.

## Value

`SplitMatchProbability()` returns a numeric giving the proportion of
permissible non-trivial splits that divide the terminals into
bipartitions of the sizes given, that match as well as `split1` and
`split2` do.

`LnSplitMatchProbability()` returns the natural logarithm of the
probability.

## See also

Other split information functions:
[`CharacterInformation()`](https://ms609.github.io/TreeTools/reference/CharacterInformation.md),
[`SplitInformation()`](https://ms609.github.io/TreeTools/reference/SplitInformation.md),
[`TreesMatchingSplit()`](https://ms609.github.io/TreeTools/reference/TreesMatchingSplit.md),
[`UnrootedTreesMatchingSplit()`](https://ms609.github.io/TreeTools/reference/UnrootedTreesMatchingSplit.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
split1 <- as.Splits(c(rep(TRUE, 4), rep(FALSE, 4)))
split2 <- as.Splits(c(rep(TRUE, 3), rep(FALSE, 5)))
SplitMatchProbability(split1, split2)
#> [1] 0.1428571
LnSplitMatchProbability(split1, split2)
#> [1] -1.94591
```
