# Double factorial

Calculate the double factorial of a number, or its logarithm.

## Usage

``` r
DoubleFactorial(n)

DoubleFactorial64(n)

LnDoubleFactorial(n)

Log2DoubleFactorial(n)

LogDoubleFactorial(n)

LnDoubleFactorial.int(n)

LogDoubleFactorial.int(n)
```

## Arguments

- n:

  Vector of integers.

## Value

Returns the double factorial, *n* \* (*n* - 2) \* (*n* - 4) \* (*n* - 6)
\* ...

## Functions

- `DoubleFactorial64()`: Returns the exact double factorial as a 64-bit
  `integer64`, for `n` \< 34.

- `LnDoubleFactorial()`: Returns the logarithm of the double factorial.

- `Log2DoubleFactorial()`: Returns the logarithm of the double
  factorial.

- `LnDoubleFactorial.int()`: Slightly faster, when x is known to be
  length one and below 50001

## See also

Other double factorials:
[`doubleFactorials`](https://ms609.github.io/TreeTools/dev/reference/doubleFactorials.md),
[`logDoubleFactorials`](https://ms609.github.io/TreeTools/dev/reference/logDoubleFactorials.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
DoubleFactorial (-4:0) # Return 1 if n < 2
#> [1] 1 1 1 1 1
DoubleFactorial (2) # 2
#> [1] 2
DoubleFactorial (5) # 1 * 3 * 5
#> [1] 15
exp(LnDoubleFactorial.int (8)) # log(2 * 4 * 6 * 8)
#> [1] 384
DoubleFactorial64(31)
#> integer64
#> [1] 191898783962510625
```
