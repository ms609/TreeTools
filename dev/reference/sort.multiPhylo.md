# Sort a list of phylogenetic trees

Trees are sorted by their [mixed base
representation](https://ms609.github.io/TreeTools/dev/reference/TreeNumber.md),
treating their leaves in the order of their labels (i.e. alphabetically,
if leaves are labelled with text).

## Usage

``` r
# S3 method for class 'multiPhylo'
sort(x, decreasing = FALSE, na.last = NA, ...)

# S3 method for class 'phylo'
e1 == e2

# S3 method for class 'phylo'
e1 < e2

# S3 method for class 'phylo'
e1 > e2

# S3 method for class 'MixedBase'
e1 == e2

# S3 method for class 'MixedBase'
e1 < e2

# S3 method for class 'MixedBase'
e1 > e2
```

## Arguments

- x, decreasing, na.last, ...:

  As in [`sort()`](https://rdrr.io/r/base/sort.html).

- e1, e2:

  Objects to be compared.

## See also

Other utility functions:
[`ClusterTable`](https://ms609.github.io/TreeTools/dev/reference/ClusterTable.md),
[`ClusterTable-methods`](https://ms609.github.io/TreeTools/dev/reference/ClusterTable-methods.md),
[`Hamming()`](https://ms609.github.io/TreeTools/dev/reference/Hamming.md),
[`MSTEdges()`](https://ms609.github.io/TreeTools/dev/reference/MSTEdges.md),
[`SampleOne()`](https://ms609.github.io/TreeTools/dev/reference/SampleOne.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/dev/reference/TipTimedTree.md),
[`UnshiftTree()`](https://ms609.github.io/TreeTools/dev/reference/UnshiftTree.md),
[`as.multiPhylo()`](https://ms609.github.io/TreeTools/dev/reference/as.multiPhylo.md),
[`match,phylo,phylo-method`](https://ms609.github.io/TreeTools/dev/reference/match.multiPhylo.md),
[`sapply64()`](https://ms609.github.io/TreeTools/dev/reference/sapply64.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
sort(as.phylo(5:0, 7))
#> 6 phylogenetic trees
```
