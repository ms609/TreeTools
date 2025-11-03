# Remove quotation marks from a string

Remove quotation marks from a string

## Usage

``` r
Unquote(string)
```

## Arguments

- string:

  Input string

## Value

`Unquote()` returns `string`, with any matched punctuation marks and
trailing whitespace removed.

## See also

Other string parsing functions:
[`EndSentence()`](https://ms609.github.io/TreeTools/reference/EndSentence.md),
[`MatchStrings()`](https://ms609.github.io/TreeTools/reference/MatchStrings.md),
[`MorphoBankDecode()`](https://ms609.github.io/TreeTools/reference/MorphoBankDecode.md),
[`RightmostCharacter()`](https://ms609.github.io/TreeTools/reference/RightmostCharacter.md)

## Author

Martin R. Smith

## Examples

``` r
Unquote("'Hello World'")
#> [1] "Hello World"
```
