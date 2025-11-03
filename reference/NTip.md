# Number of leaves in a phylogenetic tree

`NTip()` extends
[`ape::Ntip()`](https://rdrr.io/pkg/ape/man/summary.phylo.html) to
handle objects of class `Splits` and `list`, and edge matrices
(equivalent to `tree$edge`).

## Usage

``` r
NTip(phy)

# Default S3 method
NTip(phy)

# S3 method for class 'Splits'
NTip(phy)

# S3 method for class 'list'
NTip(phy)

# S3 method for class 'phylo'
NTip(phy)

# S3 method for class 'multiPhylo'
NTip(phy)

# S3 method for class 'phyDat'
NTip(phy)

# S3 method for class 'matrix'
NTip(phy)
```

## Arguments

- phy:

  Object representing one or more phylogenetic trees.

## Value

`NTip()` returns an integer specifying the number of tips in each object
in `phy`.

## See also

Other tree properties:
[`Cherries()`](https://ms609.github.io/TreeTools/reference/Cherries.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/reference/ConsensusWithout.md),
[`LongBranch()`](https://ms609.github.io/TreeTools/reference/LongBranch.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/reference/MatchEdges.md),
[`NSplits()`](https://ms609.github.io/TreeTools/reference/NSplits.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/reference/NodeNumbers.md),
[`PathLengths()`](https://ms609.github.io/TreeTools/reference/PathLengths.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/reference/TipLabels.md),
[`TreeIsRooted()`](https://ms609.github.io/TreeTools/reference/TreeIsRooted.md),
[`Treeness()`](https://ms609.github.io/TreeTools/reference/Treeness.md)

Other Splits operations:
[`LabelSplits()`](https://ms609.github.io/TreeTools/reference/LabelSplits.md),
[`NSplits()`](https://ms609.github.io/TreeTools/reference/NSplits.md),
[`PolarizeSplits()`](https://ms609.github.io/TreeTools/reference/PolarizeSplits.md),
[`SplitFrequency()`](https://ms609.github.io/TreeTools/reference/SplitFrequency.md),
[`Splits`](https://ms609.github.io/TreeTools/reference/Splits.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/reference/TipLabels.md),
[`TipsInSplits()`](https://ms609.github.io/TreeTools/reference/TipsInSplits.md),
[`match,Splits,Splits-method`](https://ms609.github.io/TreeTools/reference/match.Splits.md),
[`xor()`](https://ms609.github.io/TreeTools/reference/xor.md)
