# Split matching

`match()` returns a vector of the positions of (first) matches of splits
in its first argument in its second. `%in%` is a more intuitive
interface as a binary operator, which returns a logical vector
indicating whether there is a match or not for each split in its left
operand.

## Usage

``` r
# S4 method for class 'Splits,Splits'
match(x, table, nomatch = NA_integer_, incomparables = NULL)

match(x, table, nomatch = NA_integer_, incomparables = NULL)

# S4 method for class 'Splits,Splits'
x %in% table
```

## Arguments

- x, table:

  Object of class `Splits`.

- nomatch:

  Integer value that will be used in place of `NA` in the case where no
  match is found.

- incomparables:

  Ignored. (Included for consistency with generic.)

## Value

`match()` returns an integer vector specifying the position in `table`
that matches each element in `x`, or `nomatch` if no match is found.

## See also

Corresponding base functions are documented in
[`match()`](https://rdrr.io/r/base/match.html).

Other Splits operations:
[`LabelSplits()`](https://ms609.github.io/TreeTools/reference/LabelSplits.md),
[`NSplits()`](https://ms609.github.io/TreeTools/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/reference/NTip.md),
[`PolarizeSplits()`](https://ms609.github.io/TreeTools/reference/PolarizeSplits.md),
[`SplitFrequency()`](https://ms609.github.io/TreeTools/reference/SplitFrequency.md),
[`Splits`](https://ms609.github.io/TreeTools/reference/Splits.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/reference/TipLabels.md),
[`TipsInSplits()`](https://ms609.github.io/TreeTools/reference/TipsInSplits.md),
[`xor()`](https://ms609.github.io/TreeTools/reference/xor.md)

## Examples

``` r
splits1 <- as.Splits(BalancedTree(7))
splits2 <- as.Splits(PectinateTree(7))

match(splits1, splits2)
#> [1]  3  1 NA NA
```
