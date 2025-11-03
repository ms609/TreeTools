# Total Cophenetic Index

`TotalCopheneticIndex()` calculates the total cophenetic index (Mir et
al. 2013) for any tree, a measure of its balance; `TCIContext()` lists
its possible values.

## Usage

``` r
TotalCopheneticIndex(x)

TCIContext(x)

# S3 method for class 'numeric'
TCIContext(x)
```

## Arguments

- x:

  A tree of class `phylo`, its `$edge` property, or a list thereof.

## Value

`TotalCopheneticIndex()` returns an integer denoting the total
cophenetic index.

`TCIContext()` returns a data frame detailing the maximum and minimum
value obtainable for the Total Cophenetic Index for rooted binary trees
with the number of leaves of the given tree, and the expected value
under the Yule and Uniform models. The variance of the expected value is
given under the Yule model, but cannot be obtained by calculation for
the Uniform model.

## Details

The Total Cophenetic Index is a measure of tree balance – i.e. whether a
(phylogenetic) tree comprises symmetric pairs of nodes, or has a
pectinate "caterpillar" shape. The index has a greater resolution power
than Sackin's and Colless' indices, and can be applied to trees that are
not perfectly resolved.

For a tree with *n* leaves, the Total Cophenetic Index can take values
of 0 to `choose(n, 3)`. The minimum value is higher for a perfectly
resolved (i.e. dichotomous) tree (see Lemma 14 of Mir *et al.* 2013).
Formulae to calculate the expected values under the Yule and Uniform
models of evolution are given in Theorems 17 and 23.

Full details are provided by Mir et al. (2013) .

The J¹ index (Lemant et al. 2022) has advantages over the Total
Cophenetic Index, particularly when comparing trees with different
numbers of leaves, or where the population size of nodes is meaningful;
see
[`J1Index()`](https://ms609.github.io/TreeTools/reference/J1Index.md).

## References

Lemant J, Le Sueur C, Manojlović V, Noble R (2022). “Robust, Universal
Tree Balance Indices.” *Systematic Biology*, **71**(5), 1210–1224.
[doi:10.1093/sysbio/syac027](https://doi.org/10.1093/sysbio/syac027) .  
  
Mir A, Rosselló F, Rotger LA (2013). “A new balance index for
phylogenetic trees.” *Mathematical Biosciences*, **241**(1), 125–136.
[doi:10.1016/j.mbs.2012.10.005](https://doi.org/10.1016/j.mbs.2012.10.005)
.

## See also

- [`J1Index()`](https://ms609.github.io/TreeTools/reference/J1Index.md)
  provides a more robust, universal tree balance index.

- `cophen.index()` in the package
  [CollessLike](https://github.com/LuciaRotger/CollessLike) provides an
  alternative implementation of this index and its predecessors.

Other tree characterization functions:
[`CladisticInfo()`](https://ms609.github.io/TreeTools/reference/CladisticInfo.md),
[`Consensus()`](https://ms609.github.io/TreeTools/reference/Consensus.md),
[`J1Index()`](https://ms609.github.io/TreeTools/reference/J1Index.md),
[`Stemwardness`](https://ms609.github.io/TreeTools/reference/Stemwardness.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
# Balanced trees have the minimum index for a binary tree;
# Pectinate trees the maximum:
TCIContext(8)
#>   maximum minimum uniform.expected yule.expected yule.variance
#> 1      56      16          38.8345      28.51429      90.52281
TotalCopheneticIndex(PectinateTree(8))
#> [1] 56
TotalCopheneticIndex(BalancedTree(8))
#> [1] 16
TotalCopheneticIndex(StarTree(8))
#> [1] 0


# Examples from Mir et al. (2013):
tree12 <- ape::read.tree(text="(1, (2, (3, (4, 5))));")  #Fig. 4, tree 12
TotalCopheneticIndex(tree12) # 10
#> [1] 10
tree8  <- ape::read.tree(text="((1, 2, 3, 4), 5);")      #Fig. 4, tree 8
TotalCopheneticIndex(tree8)  # 6
#> [1] 6
TCIContext(tree8)
#>   maximum minimum uniform.expected yule.expected yule.variance
#> 1      10       5         8.285714      7.166667      5.138889
TCIContext(5L) # Context for a tree with 5 leaves.
#>   maximum minimum uniform.expected yule.expected yule.variance
#> 1      10       5         8.285714      7.166667      5.138889
```
