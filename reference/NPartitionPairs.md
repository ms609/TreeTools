# Distributions of tips consistent with a partition pair

`NPartitionPairs()` calculates the number of terminal arrangements
matching a specified configuration of two splits.

## Usage

``` r
NPartitionPairs(configuration)
```

## Arguments

- configuration:

  Integer vector of length four specifying the number of terminals that
  occur in both (1) splits A1 and A2; (2) splits A1 and B2; (3) splits
  B1 and A2; (4) splits B1 and B2.

## Value

The number of ways to distribute `sum(configuration)` taxa according to
the specified pattern.

## Details

Consider splits that divide eight terminals, labelled A to H.

|                |           |           |            |
|----------------|-----------|-----------|------------|
| Bipartition 1: | ABCD:EFGH | A1 = ABCD | B1 = EFGH  |
| Bipartition 2: | ABE:CDFGH | A2 = ABE  | B2 = CDFGH |

This can be represented by an association matrix:

|      |      |      |
|------|------|------|
|      | *A2* | *B2* |
| *A1* | AB   | C    |
| *B1* | E    | FGH  |

The cells in this matrix contain 2, 1, 1 and 3 terminals respectively;
this four-element vector (`c(2, 1, 1, 3)`) is the `configuration`
implied by this pair of bipartition splits.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
NPartitionPairs(c(2, 1, 1, 3))
#> [1] 12
```
