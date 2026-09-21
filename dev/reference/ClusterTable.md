# Convert phylogenetic tree to `ClusterTable`

`as.ClusterTable()` converts a phylogenetic tree to a `ClusterTable`
object, which is an internal representation of its splits suitable for
rapid tree distance calculation (per Day, 1985).

## Usage

``` r
as.ClusterTable(x, tipLabels = NULL, ...)

# S3 method for class 'phylo'
as.ClusterTable(x, tipLabels = NULL, ...)

# S3 method for class 'list'
as.ClusterTable(x, tipLabels = NULL, ...)

# S3 method for class 'multiPhylo'
as.ClusterTable(x, tipLabels = NULL, ...)
```

## Arguments

- x:

  Object to convert into `ClusterTable`: perhaps a tree of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- tipLabels:

  Character vector specifying sequence in which to order tip labels.

- ...:

  Unused.

## Value

`as.ClusterTable()` returns an object of class `ClusterTable`, or a list
thereof.

## Details

Each row of a cluster table relates to a clade on a tree rooted on
tip 1. Tips are numbered according to the order in which they are
visited in preorder: i.e., if plotted using `plot(x)`, from the top of
the page downwards. A clade containing the tips 2 .. 5 would be denoted
by the entry `2, 5`, in either row 2 or row 5 of the cluster table.

## References

Day WHE (1985). “Optimal algorithms for comparing trees with labeled
leaves.” *Journal of Classification*, **2**(1), 7–28.
[doi:10.1007/BF01908061](https://doi.org/10.1007/BF01908061) .

## See also

[S3
methods](https://ms609.github.io/TreeTools/dev/reference/ClusterTable-methods.md)
for `ClusterTable` objects.

Other utility functions:
[`ClusterTable-methods`](https://ms609.github.io/TreeTools/dev/reference/ClusterTable-methods.md),
[`Hamming()`](https://ms609.github.io/TreeTools/dev/reference/Hamming.md),
[`MSTEdges()`](https://ms609.github.io/TreeTools/dev/reference/MSTEdges.md),
[`SampleOne()`](https://ms609.github.io/TreeTools/dev/reference/SampleOne.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/dev/reference/TipTimedTree.md),
[`UnshiftTree()`](https://ms609.github.io/TreeTools/dev/reference/UnshiftTree.md),
[`as.multiPhylo()`](https://ms609.github.io/TreeTools/dev/reference/as.multiPhylo.md),
[`match,phylo,phylo-method`](https://ms609.github.io/TreeTools/dev/reference/match.multiPhylo.md),
[`sapply64()`](https://ms609.github.io/TreeTools/dev/reference/sapply64.md),
[`sort.multiPhylo()`](https://ms609.github.io/TreeTools/dev/reference/sort.multiPhylo.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tree1 <- ape::read.tree(text = "(A, (B, (C, (D, E))));");
tree2 <- ape::read.tree(text = "(A, (B, (D, (C, E))));");
ct1 <- as.ClusterTable(tree1)
summary(ct1)
#> ClusterTable on 5 leaves:
#>  12345
#>  **...
#>  ***..
#>  ****.
#>  *****
#>  1: E  2: D  3: C  4: B  5: A 
as.matrix(ct1)
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    1    2
#> [3,]    1    3
#> [4,]    1    4
#> [5,]    1    5

# Tip label order must match ct1 to allow comparison
ct2 <- as.ClusterTable(tree2, tipLabels = LETTERS[1:5])

# It can thus be safer to use
ctList <- as.ClusterTable(c(tree1, tree2))
ctList[[2]]
#> ClusterTable on 5 leaves: A .. E
```
