# Add full stop to end of a sentence

Add full stop to end of a sentence

## Usage

``` r
EndSentence(string)
```

## Arguments

- string:

  Input string

## Value

`EndSentence()` returns `string`, punctuated with a final full stop
(period).\`

## See also

Other string parsing functions:
[`MatchStrings()`](https://ms609.github.io/TreeTools/dev/reference/MatchStrings.md),
[`MorphoBankDecode()`](https://ms609.github.io/TreeTools/dev/reference/MorphoBankDecode.md),
[`RightmostCharacter()`](https://ms609.github.io/TreeTools/dev/reference/RightmostCharacter.md),
[`Unquote()`](https://ms609.github.io/TreeTools/dev/reference/Unquote.md)

## Author

Martin R. Smith

## Examples

``` r
EndSentence("Hello World") # "Hello World."
#> [1] "Hello World."
```
