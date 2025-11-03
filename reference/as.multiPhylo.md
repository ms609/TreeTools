# Convert object to `multiPhylo` class

Converts representations of phylogenetic trees to an object of the "ape"
class `multiPhylo`.

## Usage

``` r
as.multiPhylo(x)

# S3 method for class 'phylo'
as.multiPhylo(x)

# S3 method for class 'list'
as.multiPhylo(x)

# S3 method for class 'phyDat'
as.multiPhylo(x)

# S3 method for class 'Splits'
as.multiPhylo(x)
```

## Arguments

- x:

  Object to be converted

## Value

`as.multiPhylo` returns an object of class `multiPhylo`

`as.multiPhylo.phyDat()` returns a list of trees, each corresponding to
the partitions implied by each non-ambiguous character in `x`.

## See also

Other utility functions:
[`ClusterTable`](https://ms609.github.io/TreeTools/reference/ClusterTable.md),
[`ClusterTable-methods`](https://ms609.github.io/TreeTools/reference/ClusterTable-methods.md),
[`Hamming()`](https://ms609.github.io/TreeTools/reference/Hamming.md),
[`MSTEdges()`](https://ms609.github.io/TreeTools/reference/MSTEdges.md),
[`SampleOne()`](https://ms609.github.io/TreeTools/reference/SampleOne.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/reference/TipTimedTree.md),
[`UnshiftTree()`](https://ms609.github.io/TreeTools/reference/UnshiftTree.md),
[`match,phylo,phylo-method`](https://ms609.github.io/TreeTools/reference/match.multiPhylo.md),
[`sapply64()`](https://ms609.github.io/TreeTools/reference/sapply64.md),
[`sort.multiPhylo()`](https://ms609.github.io/TreeTools/reference/sort.multiPhylo.md)

## Examples

``` r
as.multiPhylo(BalancedTree(8))
#> 1 phylogenetic tree
as.multiPhylo(list(BalancedTree(8), PectinateTree(8)))
#> 2 phylogenetic trees
data("Lobo")
as.multiPhylo(Lobo.phy)
#> 115 phylogenetic trees
```
