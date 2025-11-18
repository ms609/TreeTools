# Identify and remove trivial splits

`TrivialSplits()` identifies trivial splits (which separate one or zero
leaves from all others); `WithoutTrivialSplits()` removes them from a
`Splits` object.

## Usage

``` r
TrivialSplits(splits, nTip = attr(splits, "nTip"))

WithoutTrivialSplits(splits, nTip = attr(splits, "nTip"))
```

## Arguments

- splits:

  An object of class
  [`Splits`](https://ms609.github.io/TreeTools/dev/reference/Splits.md).

- nTip:

  Integer specifying number of tips (leaves).

## Value

`TrivialSplits()` returns a logical vector specifying whether each split
in `splits` is trivial, i.e. includes or excludes only a single tip or
no tips at all.

`WithoutTrivialSplits()` returns a `Splits` object with trivial splits
removed.

## See also

Other split manipulation functions:
[`DropTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md),
[`SplitConsistent()`](https://ms609.github.io/TreeTools/dev/reference/SplitConsistent.md),
[`Subsplit()`](https://ms609.github.io/TreeTools/dev/reference/Subsplit.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
splits <- as.Splits(PectinateTree(letters[1:9]))
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
summary(WithoutTrivialSplits(efgh))
#> 1 bipartition split dividing 4 tips, e .. h
#>      1234
#>  16  ..**
#> 
#>  Tip 1: e     Tip 2: f    Tip 3: g    Tip 4: h   
```
