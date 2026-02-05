# Rightmost character of string

`RightmostCharacter()` is a convenience function that returns the final
character of a string.

## Usage

``` r
RightmostCharacter(string, len = nchar(string))
```

## Arguments

- string:

  Character string.

- len:

  (Optional) Integer specifying number of characters in `string`.

## Value

`RightmostCharacter()` returns the rightmost character of a string.

## See also

Other string parsing functions:
[`EndSentence()`](https://ms609.github.io/TreeTools/reference/EndSentence.md),
[`MatchStrings()`](https://ms609.github.io/TreeTools/reference/MatchStrings.md),
[`MorphoBankDecode()`](https://ms609.github.io/TreeTools/reference/MorphoBankDecode.md),
[`Unquote()`](https://ms609.github.io/TreeTools/reference/Unquote.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
RightmostCharacter("Hello, World!")
#> [1] "!"
```
