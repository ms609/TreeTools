# Check for mismatch between character vectors

Checks that entries in one character vector occur in another, suggesting
corrections for mismatched elements.

## Usage

``` r
MatchStrings(x, table, Fail = stop, max.distance = 0.5, ...)
```

## Arguments

- x, table:

  Character vectors, in which all elements of `x` are expected to occur
  in `table`.

- Fail:

  Function to call if a mismatch is found.

- max.distance, ...:

  Arguments to [`agrep()`](https://rdrr.io/r/base/agrep.html), used to
  propose possible matches to the user.

## Value

`MatchStrings()` returns the elements of `x` that occur in `table`.

## See also

Other string parsing functions:
[`EndSentence()`](https://ms609.github.io/TreeTools/dev/reference/EndSentence.md),
[`MorphoBankDecode()`](https://ms609.github.io/TreeTools/dev/reference/MorphoBankDecode.md),
[`RightmostCharacter()`](https://ms609.github.io/TreeTools/dev/reference/RightmostCharacter.md),
[`Unquote()`](https://ms609.github.io/TreeTools/dev/reference/Unquote.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tree <- BalancedTree(8)
MatchStrings(c("t1", "tip2", "t3"), TipLabels(tree), Fail = message)
#> Could not find 'tip2' in TipLabels(tree).  Did you mean 't2'?
#> [1] "t1" "t3"
```
