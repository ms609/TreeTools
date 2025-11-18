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
[`Cherries()`](https://ms609.github.io/TreeTools/dev/reference/Cherries.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/dev/reference/ConsensusWithout.md),
[`LongBranch()`](https://ms609.github.io/TreeTools/dev/reference/LongBranch.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/dev/reference/MatchEdges.md),
[`NSplits()`](https://ms609.github.io/TreeTools/dev/reference/NSplits.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/dev/reference/NodeNumbers.md),
[`PathLengths()`](https://ms609.github.io/TreeTools/dev/reference/PathLengths.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/dev/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/dev/reference/TipLabels.md),
[`TreeIsRooted()`](https://ms609.github.io/TreeTools/dev/reference/TreeIsRooted.md),
[`Treeness()`](https://ms609.github.io/TreeTools/dev/reference/Treeness.md)

Other Splits operations:
[`LabelSplits()`](https://ms609.github.io/TreeTools/dev/reference/LabelSplits.md),
[`NSplits()`](https://ms609.github.io/TreeTools/dev/reference/NSplits.md),
[`PolarizeSplits()`](https://ms609.github.io/TreeTools/dev/reference/PolarizeSplits.md),
[`SplitFrequency()`](https://ms609.github.io/TreeTools/dev/reference/SplitFrequency.md),
[`Splits`](https://ms609.github.io/TreeTools/dev/reference/Splits.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/dev/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/dev/reference/TipLabels.md),
[`TipsInSplits()`](https://ms609.github.io/TreeTools/dev/reference/TipsInSplits.md),
[`match,Splits,Splits-method`](https://ms609.github.io/TreeTools/dev/reference/match.Splits.md),
[`xor()`](https://ms609.github.io/TreeTools/dev/reference/xor.md)
