# Number of trees one SPR step away

`N1Spr()` calculates the number of trees one subtree prune-and-regraft
operation away from a binary input tree using the formula given by Allen
and Steel (2001) ; `IC1Spr()` calculates the information content of
trees at this distance: i.e. the entropy corresponding to the proportion
of all possible *n*-tip trees whose SPR distance is at most one from a
specified tree.

## Usage

``` r
N1Spr(n)

IC1Spr(n)
```

## Arguments

- n:

  Integer vector specifying the number of tips in a tree.

## Value

`N1Spr()` returns an integer vector denoting the number of trees one SPR
rearrangement away from the input tree..

`IC1Spr()` returns an numeric vector giving the phylogenetic information
content of trees 0 or 1 SPR rearrangement from an *n*-leaf tree, in
bits.

## References

Allen BL, Steel MA (2001). “Subtree transfer operations and their
induced metrics on evolutionary trees.” *Annals of Combinatorics*,
**5**(1), 1–15.
[doi:10.1007/s00026-001-8006-8](https://doi.org/10.1007/s00026-001-8006-8)
.

## Examples

``` r
N1Spr(4:6)
#> [1]  2 12 30
IC1Spr(5)
#> [1] 0.2064509
```
