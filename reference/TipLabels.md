# Extract tip labels

`TipLabels()` extracts labels from an object: for example, names of taxa
in a phylogenetic tree or data matrix. `AllTipLabels()` extracts all
labels, where entries of a list of trees may pertain to different taxa.

## Usage

``` r
TipLabels(x, single = TRUE)

# Default S3 method
TipLabels(x, single = TRUE)

# S3 method for class 'matrix'
TipLabels(x, single = TRUE)

# S3 method for class 'logical'
TipLabels(x, single = TRUE)

# S3 method for class 'phylo'
TipLabels(x, single = TRUE)

# S3 method for class 'phyDat'
TipLabels(x, single = TRUE)

# S3 method for class 'MixedBase'
TipLabels(x, single = TRUE)

# S3 method for class 'TreeNumber'
TipLabels(x, single = TRUE)

# S3 method for class 'Splits'
TipLabels(x, single = TRUE)

# S3 method for class 'list'
TipLabels(x, single = FALSE)

# S3 method for class 'multiPhylo'
TipLabels(x, single = FALSE)

# S3 method for class 'character'
TipLabels(x, single = TRUE)

# S3 method for class 'numeric'
TipLabels(x, single = TRUE)

# S3 method for class 'phyDat'
TipLabels(x, single = TRUE)

AllTipLabels(x)

# S3 method for class 'list'
AllTipLabels(x)

# S3 method for class 'multiPhylo'
AllTipLabels(x)

# S3 method for class 'phylo'
AllTipLabels(x)

# S3 method for class 'Splits'
AllTipLabels(x)

# S3 method for class 'TreeNumber'
AllTipLabels(x)

# S3 method for class 'matrix'
AllTipLabels(x)
```

## Arguments

- x:

  An object of a supported class (see Usage section above).

- single:

  Logical specifying whether to report the labels for the first object
  only (`TRUE`), or for each object in a list (`FALSE`).

## Value

`TipLabels()` returns a character vector listing the tip labels
appropriate to `x`. If `x` is a single integer, this will be a vector
`t1`, `t2` ... `tx`, to match the default of
[`rtree()`](https://rdrr.io/pkg/ape/man/rtree.html).

## See also

Other tree properties:
[`Cherries()`](https://ms609.github.io/TreeTools/reference/Cherries.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/reference/ConsensusWithout.md),
[`LongBranch()`](https://ms609.github.io/TreeTools/reference/LongBranch.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/reference/MatchEdges.md),
[`NSplits()`](https://ms609.github.io/TreeTools/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/reference/NTip.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/reference/NodeNumbers.md),
[`PathLengths()`](https://ms609.github.io/TreeTools/reference/PathLengths.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/reference/SplitsInBinaryTree.md),
[`TreeIsRooted()`](https://ms609.github.io/TreeTools/reference/TreeIsRooted.md),
[`Treeness()`](https://ms609.github.io/TreeTools/reference/Treeness.md)

Other Splits operations:
[`LabelSplits()`](https://ms609.github.io/TreeTools/reference/LabelSplits.md),
[`NSplits()`](https://ms609.github.io/TreeTools/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/reference/NTip.md),
[`PolarizeSplits()`](https://ms609.github.io/TreeTools/reference/PolarizeSplits.md),
[`SplitFrequency()`](https://ms609.github.io/TreeTools/reference/SplitFrequency.md),
[`Splits`](https://ms609.github.io/TreeTools/reference/Splits.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/reference/SplitsInBinaryTree.md),
[`TipsInSplits()`](https://ms609.github.io/TreeTools/reference/TipsInSplits.md),
[`match,Splits,Splits-method`](https://ms609.github.io/TreeTools/reference/match.Splits.md),
[`xor()`](https://ms609.github.io/TreeTools/reference/xor.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
TipLabels(BalancedTree(letters[5:1]))
#> [1] "e" "d" "c" "b" "a"
TipLabels(5)
#> [1] "t1" "t2" "t3" "t4" "t5"

data("Lobo")
head(TipLabels(Lobo.phy))
#> [1] "Tubiluchus_Priapulida"  "Cricocosmia"            "Aysheaia"              
#> [4] "Siberion"               "Onychodictyon_ferox"    "Onychodictyon_gracilis"

AllTipLabels(c(BalancedTree(4), PectinateTree(8)))
#> [1] "t1" "t2" "t3" "t4" "t5" "t6" "t7" "t8"
```
