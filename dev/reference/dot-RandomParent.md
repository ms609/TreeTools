# Random parent vector

Random parent vector

## Usage

``` r
.RandomParent(n, seed = sample.int(2147483647L, 1L))
```

## Arguments

- n:

  Integer specifying number of leaves.

- seed:

  (Optional) Integer with which to seed Mersenne Twister random number
  generator in C++.

## Value

Integer vector corresponding to the "parent" entry of `tree[["edge"]]`,
where the "child" entry, i.e. column 2, is numbered sequentially from
`1:n`.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)
