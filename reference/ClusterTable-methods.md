# S3 methods for `ClusterTable` objects

S3 methods for
[`ClusterTable`](https://ms609.github.io/TreeTools/reference/ClusterTable.md)
objects.

## Usage

``` r
# S3 method for class 'ClusterTable'
as.matrix(x, ...)

# S3 method for class 'ClusterTable'
print(x, ...)

# S3 method for class 'ClusterTable'
summary(object, ...)
```

## Arguments

- x, object:

  Object of class `ClusterTable`.

- ...:

  Additional arguments for consistency with S3 methods.

## See also

Other utility functions:
[`ClusterTable`](https://ms609.github.io/TreeTools/reference/ClusterTable.md),
[`Hamming()`](https://ms609.github.io/TreeTools/reference/Hamming.md),
[`MSTEdges()`](https://ms609.github.io/TreeTools/reference/MSTEdges.md),
[`SampleOne()`](https://ms609.github.io/TreeTools/reference/SampleOne.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/reference/TipTimedTree.md),
[`UnshiftTree()`](https://ms609.github.io/TreeTools/reference/UnshiftTree.md),
[`as.multiPhylo()`](https://ms609.github.io/TreeTools/reference/as.multiPhylo.md),
[`match,phylo,phylo-method`](https://ms609.github.io/TreeTools/reference/match.multiPhylo.md),
[`sapply64()`](https://ms609.github.io/TreeTools/reference/sapply64.md),
[`sort.multiPhylo()`](https://ms609.github.io/TreeTools/reference/sort.multiPhylo.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
clustab <- as.ClusterTable(TreeTools::BalancedTree(6))
as.matrix(clustab)
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    2    3
#> [3,]    1    3
#> [4,]    1    4
#> [5,]    1    5
#> [6,]    1    6
print(clustab)
#> ClusterTable on 6 leaves: t1 .. t6
summary(clustab)
#> ClusterTable on 6 leaves:
#>  123456
#>  .**...
#>  ***...
#>  ****..
#>  *****.
#>  ******
#>  1: t6  2: t5  3: t4  4: t3  5: t2  6: t1 
```
