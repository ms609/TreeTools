# Apply a function that returns 64-bit integers over a list or vector

Wrappers for members of the
[`lapply()`](https://rdrr.io/r/base/lapply.html) family intended for use
when a function `FUN` returns a vector of `integer64` objects.
[`vapply()`](https://rdrr.io/r/base/lapply.html),
[`sapply()`](https://rdrr.io/r/base/lapply.html) or
[`replicate()`](https://rdrr.io/r/base/lapply.html) drop the `integer64`
class, resulting in a vector of numerics that require conversion back to
64-bit integers. These functions restore the missing `class` attribute.

## Usage

``` r
sapply64(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE)

vapply64(X, FUN, FUN.LEN = 1, ...)

replicate64(n, expr, simplify = "array")
```

## Arguments

- X:

  a vector (atomic or list) or an
  [`expression`](https://rdrr.io/r/base/expression.html) object. Other
  objects (including classed objects) will be coerced by
  `base::`[`as.list`](https://rdrr.io/r/base/list.html).

- FUN:

  the function to be applied to each element of `X`: see ‘Details’. In
  the case of functions like `+`, `%*%`, the function name must be
  backquoted or quoted.

- ...:

  optional arguments to `FUN`.

- simplify:

  logical or character string; should the result be simplified to a
  vector, matrix or higher dimensional array if possible? For `sapply`
  it must be named and not abbreviated. The default value, `TRUE`,
  returns a vector or matrix if appropriate, whereas if
  `simplify = "array"` the result may be an
  [`array`](https://rdrr.io/r/base/array.html) of “rank”
  (\\=\\`length(dim(.))`) one higher than the result of `FUN(X[[i]])`.

- USE.NAMES:

  logical; if `TRUE` and if `X` is character, use `X` as
  [`names`](https://rdrr.io/r/base/names.html) for the result unless it
  had names already. Since this argument follows `...` its name cannot
  be abbreviated.

- FUN.LEN:

  Integer specifying the length of the output of `FUN`.

- n:

  integer: the number of replications.

- expr:

  the expression (a [language
  object](https://rdrr.io/r/base/is.language.html), usually a call) to
  evaluate repeatedly.

## Details

For details of the underlying functions, see
[`base::lapply()`](https://rdrr.io/r/base/lapply.html).

## See also

[`integer64()`](https://rdrr.io/pkg/bit64/man/bit64-package.html)

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
[`sort.multiPhylo()`](https://ms609.github.io/TreeTools/dev/reference/sort.multiPhylo.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
sapply64(as.phylo(1:6, 6), as.TreeNumber)
#> integer64
#> [1] 1 2 3 4 5 6
vapply64(as.phylo(1:6, 6), as.TreeNumber, 1)
#> integer64
#> [1] 1 2 3 4 5 6
set.seed(0)
replicate64(6, as.TreeNumber(RandomTree(6)))
#> integer64
#> [1] 91  45  102 24  92  48 
```
