# Is an object a `TreeNumber` object?

Is an object a `TreeNumber` object?

## Usage

``` r
is.TreeNumber(x)
```

## Arguments

- x:

  R object.

## Value

`is.TreeNumber()` returns a logical vector of length one specifying
whether `x` inherits the class `"TreeNumber"`.

## See also

Other 'TreeNumber' utilities:
[`TreeNumber`](https://ms609.github.io/TreeTools/reference/TreeNumber.md),
[`print.TreeNumber()`](https://ms609.github.io/TreeTools/reference/print.TreeNumber.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
is.TreeNumber(FALSE) # FALSE 
#> [1] FALSE
is.TreeNumber(as.TreeNumber(BalancedTree(5))) # TRUE
#> [1] TRUE
```
