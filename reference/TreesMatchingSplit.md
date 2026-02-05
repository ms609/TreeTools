# Number of trees matching a bipartition split

Calculates the number of unrooted bifurcated trees that are consistent
with a bipartition split that divides taxa into groups of size `A` and
`B`.

## Usage

``` r
TreesMatchingSplit(A, B = A[2])

LnTreesMatchingSplit(A, B = A[2])

Log2TreesMatchingSplit(A, B = A[2])
```

## Arguments

- A, B:

  Integer specifying the number of taxa in each partition.

## Value

`TreesMatchingSplit()` returns a numeric specifying the number of trees
that are compatible with the given split.

`LnTreesMatchingSplit()` and `Log2TreesMatchingSplit()` give the natural
and base-2 logarithms of this number.

## See also

Other split information functions:
[`CharacterInformation()`](https://ms609.github.io/TreeTools/reference/CharacterInformation.md),
[`SplitInformation()`](https://ms609.github.io/TreeTools/reference/SplitInformation.md),
[`SplitMatchProbability()`](https://ms609.github.io/TreeTools/reference/SplitMatchProbability.md),
[`UnrootedTreesMatchingSplit()`](https://ms609.github.io/TreeTools/reference/UnrootedTreesMatchingSplit.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
TreesMatchingSplit(5, 6)
#> [1] 99225
LnTreesMatchingSplit(5, 6)
#> [1] 11.50515
Log2TreesMatchingSplit(5, 6)
#> [1] 16.59842
```
