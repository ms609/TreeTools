# Tips contained within splits

`TipsInSplits()` specifies the number of tips that occur within each
bipartition split in a `Splits` object.

## Usage

``` r
TipsInSplits(splits, keep.names = TRUE, smallest = FALSE, ...)

# S3 method for class 'Splits'
TipsInSplits(splits, keep.names = TRUE, smallest = FALSE, ...)

# S3 method for class 'phylo'
TipsInSplits(splits, keep.names = TRUE, smallest = FALSE, ...)

SplitImbalance(splits, keep.names = TRUE, ...)

# S3 method for class 'Splits'
SplitImbalance(splits, keep.names = TRUE, ...)

# S3 method for class 'phylo'
SplitImbalance(splits, keep.names = TRUE, ...)
```

## Arguments

- splits:

  Object of class `Splits` or `phylo`.

- keep.names:

  Logical specifying whether to include the names of `splits` in the
  output.

- smallest:

  Logical; if `TRUE`, return the number of leaves in the smaller
  bipartition.

- ...:

  Additional parameters to pass to
  [`as.Splits()`](https://ms609.github.io/TreeTools/dev/reference/Splits.md).

## Value

`TipsInSplits()` returns a named vector of integers, specifying the
number of tips contained within each split in `splits`.

`SplitImbalance()` returns a named vector of integers, specifying the
number of leaves within a split that are not "balanced" by a leaf
outside it; i.e. a split that divides leaves evenly has an imbalance of
zero; one that splits two tips from ten has an imbalance of 10 - 2 = 8.

## See also

Other Splits operations:
[`LabelSplits()`](https://ms609.github.io/TreeTools/dev/reference/LabelSplits.md),
[`NSplits()`](https://ms609.github.io/TreeTools/dev/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/dev/reference/NTip.md),
[`PolarizeSplits()`](https://ms609.github.io/TreeTools/dev/reference/PolarizeSplits.md),
[`SplitFrequency()`](https://ms609.github.io/TreeTools/dev/reference/SplitFrequency.md),
[`Splits`](https://ms609.github.io/TreeTools/dev/reference/Splits.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/dev/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/dev/reference/TipLabels.md),
[`match,Splits,Splits-method`](https://ms609.github.io/TreeTools/dev/reference/match.Splits.md),
[`xor()`](https://ms609.github.io/TreeTools/dev/reference/xor.md)

## Examples

``` r
tree <- PectinateTree(8)
splits <- as.Splits(tree)
TipsInSplits(splits)
#> 11 12 13 14 15 
#>  6  5  4  3  2 

plot(tree)
LabelSplits(tree, as.character(splits), frame = "none", pos = 3L, cex = 0.7)
LabelSplits(tree, TipsInSplits(splits), unit = " tips", frame = "none",
            pos = 1L)

```
