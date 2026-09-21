# Select element at random

`SampleOne()` is a fast alternative to
[`sample()`](https://rdrr.io/r/base/sample.html) that avoids some
checks.

## Usage

``` r
SampleOne(x, len = length(x))
```

## Arguments

- x:

  A vector to sample.

- len:

  (Optional) Integer specifying length of `x`.

## Value

`SampleOne()` returns a length one vector, randomly sampled from `x`.

## See also

Other utility functions:
[`ClusterTable`](https://ms609.github.io/TreeTools/dev/reference/ClusterTable.md),
[`ClusterTable-methods`](https://ms609.github.io/TreeTools/dev/reference/ClusterTable-methods.md),
[`Hamming()`](https://ms609.github.io/TreeTools/dev/reference/Hamming.md),
[`MSTEdges()`](https://ms609.github.io/TreeTools/dev/reference/MSTEdges.md),
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
SampleOne(9:10)
#> [1] 9
SampleOne(letters[1:4])
#> [1] "a"
```
