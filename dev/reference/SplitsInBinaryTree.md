# Maximum splits in an *n*-leaf tree

`SplitsInBinaryTree()` is a convenience function to calculate the number
of splits in a fully-resolved (binary) tree with *n* leaves.

## Usage

``` r
SplitsInBinaryTree(tree)

# S3 method for class 'list'
SplitsInBinaryTree(tree)

# S3 method for class 'multiPhylo'
SplitsInBinaryTree(tree)

# S3 method for class 'numeric'
SplitsInBinaryTree(tree)

# S3 method for class '`NULL`'
SplitsInBinaryTree(tree)

# Default S3 method
SplitsInBinaryTree(tree)

# S3 method for class 'Splits'
SplitsInBinaryTree(tree)

# S3 method for class 'phylo'
SplitsInBinaryTree(tree)
```

## Arguments

- tree:

  An object of a supported format that represents a tree or set of
  trees, from which the number of leaves will be calculated.

## Value

`SplitsInBinaryTree()` returns an integer vector detailing the number of
unique non-trivial splits in a binary tree with *n* leaves.

## See also

Other tree properties:
[`Cherries()`](https://ms609.github.io/TreeTools/dev/reference/Cherries.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/dev/reference/ConsensusWithout.md),
[`LongBranch()`](https://ms609.github.io/TreeTools/dev/reference/LongBranch.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/dev/reference/MatchEdges.md),
[`NSplits()`](https://ms609.github.io/TreeTools/dev/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/dev/reference/NTip.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/dev/reference/NodeNumbers.md),
[`PathLengths()`](https://ms609.github.io/TreeTools/dev/reference/PathLengths.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/dev/reference/TipLabels.md),
[`TreeIsRooted()`](https://ms609.github.io/TreeTools/dev/reference/TreeIsRooted.md),
[`Treeness()`](https://ms609.github.io/TreeTools/dev/reference/Treeness.md)

Other Splits operations:
[`LabelSplits()`](https://ms609.github.io/TreeTools/dev/reference/LabelSplits.md),
[`NSplits()`](https://ms609.github.io/TreeTools/dev/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/dev/reference/NTip.md),
[`PolarizeSplits()`](https://ms609.github.io/TreeTools/dev/reference/PolarizeSplits.md),
[`SplitFrequency()`](https://ms609.github.io/TreeTools/dev/reference/SplitFrequency.md),
[`Splits`](https://ms609.github.io/TreeTools/dev/reference/Splits.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/dev/reference/TipLabels.md),
[`TipsInSplits()`](https://ms609.github.io/TreeTools/dev/reference/TipsInSplits.md),
[`match,Splits,Splits-method`](https://ms609.github.io/TreeTools/dev/reference/match.Splits.md),
[`xor()`](https://ms609.github.io/TreeTools/dev/reference/xor.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tree <- BalancedTree(8)
SplitsInBinaryTree(tree)
#> [1] 5
SplitsInBinaryTree(as.Splits(tree))
#> [1] 5
SplitsInBinaryTree(8)
#> [1] 5
SplitsInBinaryTree(list(tree, tree))
#> [1] 5 5
```
