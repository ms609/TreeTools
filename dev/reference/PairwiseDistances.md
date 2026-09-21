# Distances between each pair of trees

Distances between each pair of trees

## Usage

``` r
PairwiseDistances(trees, Func, valueLength = 1L, ...)
```

## Arguments

- trees:

  List of trees of class `phylo`.

- Func:

  Function returning a distance between two trees.

- valueLength:

  Integer specifying expected length of the value returned by `Func`.

- ...:

  Additional arguments to `Func`.

## Value

Matrix detailing distance between each pair of trees. Identical trees
are assumed to have zero distance.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
trees <- list(BalancedTree(8), PectinateTree(8), StarTree(8))
TCIDiff <- function(tree1, tree2) {
  TotalCopheneticIndex(tree1) - TotalCopheneticIndex(tree2)
}
PairwiseDistances(trees, TCIDiff, 1)
#>     1   2   3
#> 1     -40  16
#> 2 -40      56
#> 3  16  56    
TCIRange <- function(tree1, tree2) {
  range(TotalCopheneticIndex(tree1), TotalCopheneticIndex(tree2))
}
PairwiseDistances(trees, TCIRange, 2)
#> [[1]]
#>    1  2  3
#> 1    16  0
#> 2 16     0
#> 3  0  0   
#> 
#> [[2]]
#>    1  2  3
#> 1    56 16
#> 2 56    56
#> 3 16 56   
#> 
```
