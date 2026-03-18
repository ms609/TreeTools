# Count cherries in a tree

`Cherries()` counts the number of vertices in a binary tree whose
children are both leaves.

## Usage

``` r
Cherries(tree, nTip)

# S3 method for class 'phylo'
Cherries(tree, nTip = NTip(tree))

# S3 method for class 'numeric'
Cherries(tree, nTip)
```

## Arguments

- tree:

  A binary tree, of class `phylo`; or a matrix corresponding to its edge
  matrix.

- nTip:

  Number of leaves in tree.

## Value

`Cherries()` returns an integer specifying the number of nodes whose
children are both leaves.

Equations for the number of cherries expected under uniform and Yule
tree models are derived by McKenzie and Steel (2000) .

## References

McKenzie A, Steel M (2000). “Distributions of Cherries for Two Models of
Trees.” *Mathematical Biosciences*, **164**(1), 81–92.
[doi:10.1016/S0025-5564(99)00060-7](https://doi.org/10.1016/S0025-5564%2899%2900060-7)
.

## See also

Other tree properties:
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/reference/ConsensusWithout.md),
[`EdgeRatio()`](https://ms609.github.io/TreeTools/reference/EdgeRatio.md),
[`LongBranch()`](https://ms609.github.io/TreeTools/reference/LongBranch.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/reference/MatchEdges.md),
[`NSplits()`](https://ms609.github.io/TreeTools/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/reference/NTip.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/reference/NodeNumbers.md),
[`PathLengths()`](https://ms609.github.io/TreeTools/reference/PathLengths.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/reference/TipLabels.md),
[`TreeIsRooted()`](https://ms609.github.io/TreeTools/reference/TreeIsRooted.md),
[`Treeness()`](https://ms609.github.io/TreeTools/reference/Treeness.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)
