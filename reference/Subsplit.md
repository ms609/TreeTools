# Subset of a split on fewer leaves

`Subsplit()` removes leaves from a `Splits` object.

## Usage

``` r
Subsplit(splits, tips, keepAll = FALSE, unique = TRUE)
```

## Arguments

- splits:

  An object of class
  [`Splits`](https://ms609.github.io/TreeTools/reference/Splits.md).

- tips:

  A vector specifying a subset of the leaf labels applied to `split`.

- keepAll:

  logical specifying whether to keep entries that define trivial splits
  (i.e. splits of zero or one leaf) on the subset of leaves.

- unique:

  logical specifying whether to remove duplicate splits.

## Value

`Subsplit()` returns an object of class `Splits`, defined on the leaves
`tips`.

## See also

[`KeepTip()`](https://ms609.github.io/TreeTools/reference/DropTip.md) is
a less flexible but faster equivalent.

Other split manipulation functions:
[`DropTip()`](https://ms609.github.io/TreeTools/reference/DropTip.md),
[`SplitConsistent()`](https://ms609.github.io/TreeTools/reference/SplitConsistent.md),
[`TrivialSplits()`](https://ms609.github.io/TreeTools/reference/TrivialSplits.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
splits <- as.Splits(PectinateTree(letters[1:9]))
splits
#> 6 bipartition splits dividing 9 tips, a .. i
efgh <- Subsplit(splits, tips = letters[5:8], keepAll = TRUE)
summary(efgh)
#> 4 bipartition splits (3 trivial) dividing 4 tips, e .. h
#>      1234
#>  12  ****
#>  15  .***
#>  16  ..**
#>  17  ...*
#> 
#>  Tip 1: e     Tip 2: f    Tip 3: g    Tip 4: h   

TrivialSplits(efgh)
#>    12    15    16    17 
#>  TRUE  TRUE FALSE  TRUE 

summary(Subsplit(splits, tips = letters[5:8], keepAll = FALSE))
#> 1 bipartition split dividing 4 tips, e .. h
#>      1234
#>  16  ..**
#> 
#>  Tip 1: e     Tip 2: f    Tip 3: g    Tip 4: h   
```
