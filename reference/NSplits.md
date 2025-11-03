# Number of distinct splits

`NSplits()` counts the unique bipartition splits in a tree or object.

## Usage

``` r
NSplits(x)

NPartitions(x)

# S3 method for class 'phylo'
NSplits(x)

# S3 method for class 'list'
NSplits(x)

# S3 method for class 'multiPhylo'
NSplits(x)

# S3 method for class 'Splits'
NSplits(x)

# S3 method for class 'numeric'
NSplits(x)

# S3 method for class '`NULL`'
NSplits(x)

# S3 method for class 'ClusterTable'
NSplits(x)

# S3 method for class 'character'
NSplits(x)
```

## Arguments

- x:

  A phylogenetic tree of class `phylo`; a list of such trees (of class
  `list` or `multiPhylo`); a `Splits` object; a vector of integers; or a
  character vector listing tips of a tree, or a character of length one
  specifying a tree in Newick format.

## Value

`NSplits()` returns an integer specifying the number of bipartitions in
the specified objects, or in a binary tree with `x` tips.

## See also

Other tree properties:
[`Cherries()`](https://ms609.github.io/TreeTools/reference/Cherries.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/reference/ConsensusWithout.md),
[`LongBranch()`](https://ms609.github.io/TreeTools/reference/LongBranch.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/reference/MatchEdges.md),
[`NTip()`](https://ms609.github.io/TreeTools/reference/NTip.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/reference/NodeNumbers.md),
[`PathLengths()`](https://ms609.github.io/TreeTools/reference/PathLengths.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/reference/TipLabels.md),
[`TreeIsRooted()`](https://ms609.github.io/TreeTools/reference/TreeIsRooted.md),
[`Treeness()`](https://ms609.github.io/TreeTools/reference/Treeness.md)

Other Splits operations:
[`LabelSplits()`](https://ms609.github.io/TreeTools/reference/LabelSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/reference/NTip.md),
[`PolarizeSplits()`](https://ms609.github.io/TreeTools/reference/PolarizeSplits.md),
[`SplitFrequency()`](https://ms609.github.io/TreeTools/reference/SplitFrequency.md),
[`Splits`](https://ms609.github.io/TreeTools/reference/Splits.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/reference/TipLabels.md),
[`TipsInSplits()`](https://ms609.github.io/TreeTools/reference/TipsInSplits.md),
[`match,Splits,Splits-method`](https://ms609.github.io/TreeTools/reference/match.Splits.md),
[`xor()`](https://ms609.github.io/TreeTools/reference/xor.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
NSplits(8L)
#> [1] 5
NSplits(PectinateTree(8))
#> [1] 5
NSplits(as.Splits(BalancedTree(8)))
#> [1] 5
```
