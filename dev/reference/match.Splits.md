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

FirstMatchingSplit(x, table, nomatch, return = c("x", "table", "both"))
```

## Arguments

- x, table:

  Splits objects

- nomatch:

  Integer value that will be used in place of `NA` in the case where no
  match is found.

- incomparables:

  Ignored. (Included for consistency with generic.)

- return:

  Which index to return: in `x`, in `table`, or both

## Value

`match()` returns an integer vector specifying the position in `table`
that matches each element in `x`, or `nomatch` if no match is found.

`FirstMatchingSplit()` returns an integer (or length-2 integer if
`return = "both"`) specifying the first split in `x` to have a match in
`table` (`return = "x"`), or the index of that match
(`return = "table"`). `nomatch` (default `0`) is returned in the absence
of a match.

## See also

Corresponding base functions are documented in
[`match()`](https://rdrr.io/r/base/match.html).

Other Splits operations:
[`LabelSplits()`](https://ms609.github.io/TreeTools/dev/reference/LabelSplits.md),
[`NSplits()`](https://ms609.github.io/TreeTools/dev/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/dev/reference/NTip.md),
[`PolarizeSplits()`](https://ms609.github.io/TreeTools/dev/reference/PolarizeSplits.md),
[`SplitFrequency()`](https://ms609.github.io/TreeTools/dev/reference/SplitFrequency.md),
[`Splits`](https://ms609.github.io/TreeTools/dev/reference/Splits.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/dev/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/dev/reference/TipLabels.md),
[`TipsInSplits()`](https://ms609.github.io/TreeTools/dev/reference/TipsInSplits.md),
[`xor()`](https://ms609.github.io/TreeTools/dev/reference/xor.md)

## Examples

``` r
splits1 <- as.Splits(BalancedTree(7))
splits2 <- as.Splits(PectinateTree(7))

match(splits1, splits2)
#> [1]  3  1 NA NA
```
