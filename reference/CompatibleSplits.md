# Which splits are compatible?

Which splits are compatible?

## Usage

``` r
CompatibleSplits(splits, splits2)

.CompatibleSplit(a, b, nTip)

.CompatibleRaws(rawA, rawB, bitmask)
```

## Arguments

- splits:

  An object of class
  [`Splits`](https://ms609.github.io/TreeTools/reference/Splits.md).

- splits2:

  A second `Splits` object.

- a, b:

  [Raw](https://rdrr.io/r/base/raw.html) representations of splits, from
  a row of a `Splits` object.

- rawA, rawB:

  Raw representations of splits.

- bitmask:

  Raw masking bits that do not correspond to tips.

## Value

`CompatibleSplits` returns a logical matrix specifying whether each
split in `splits` is compatible with each split in `splits2`.

`.CompatibleSplit` returns a logical vector stating whether splits are
compatible.

`.CompatibleRaws` returns a logical vector specifying whether input raws
are compatible.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
splits <- as.Splits(BalancedTree(8))
splits2 <- as.Splits(PectinateTree(8))

summary(splits)
#> 5 bipartition splits dividing 8 tips, t1 .. t8
#>      12345678
#>  10  ****....
#>  11  **......
#>  12  ..**....
#>  14  ....**..
#>  15  ......**
#> 
#>  Tip 1: t1    Tip 2: t2   Tip 3: t3   Tip 4: t4   Tip 5: t5  
#>  Tip 6: t6    Tip 7: t7   Tip 8: t8  
summary(splits2)
#> 5 bipartition splits dividing 8 tips, t1 .. t8
#>      12345678
#>  11  ..******
#>  12  ...*****
#>  13  ....****
#>  14  .....***
#>  15  ......**
#> 
#>  Tip 1: t1    Tip 2: t2   Tip 3: t3   Tip 4: t4   Tip 5: t5  
#>  Tip 6: t6    Tip 7: t7   Tip 8: t8  

CompatibleSplits(splits, splits2)
#>      11    12   13    14   15
#> 10 TRUE  TRUE TRUE  TRUE TRUE
#> 11 TRUE  TRUE TRUE  TRUE TRUE
#> 12 TRUE FALSE TRUE  TRUE TRUE
#> 14 TRUE  TRUE TRUE FALSE TRUE
#> 15 TRUE  TRUE TRUE  TRUE TRUE
```
