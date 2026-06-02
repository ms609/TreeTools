# Convert Nexus token matrix to integer

`NexusTokensToInteger()` converts the character matrix returned by
[`ReadCharacters()`](https://ms609.github.io/TreeTools/reference/ReadCharacters.md)
to an integer matrix, mapping polymorphic, ambiguous (`?`), and
inapplicable (`-`) tokens to `NA_integer_` or to the first/last state
listed in the polymorphism, depending on `polymorphism`.

## Usage

``` r
NexusTokensToInteger(tokens, polymorphism = c("?", "first", "last"))
```

## Arguments

- tokens:

  Character matrix as returned by
  [`ReadCharacters()`](https://ms609.github.io/TreeTools/reference/ReadCharacters.md),
  a character vector as returned by
  [`NexusTokens()`](https://ms609.github.io/TreeTools/reference/ExtractTaxa.md),
  or a `phyDat` object.

- polymorphism:

  Character string specifying how to handle polymorphic tokens such as
  `"(01)"` or `"{12}"`:

  `"?"` (default)

  :   Treat as the NEXUS missing-data token: map to `NA_integer_`.

  `"first"`

  :   Use the first state digit inside the brackets.

  `"last"`

  :   Use the last state digit inside the brackets.

  Tokens `"?"` and `"-"` always map to `NA_integer_` regardless of
  `polymorphism`.

## Value

An integer matrix (or vector) with the same dimensions and `dimnames` as
`tokens`.

## Details

Only digit states `0`..`9` are recognised; non-digit symbols (and any
token whose interior contains no digits) become `NA_integer_`.
Polymorphism extraction (`polymorphism = "first"`/`"last"`) likewise
considers digits only.

If `tokens` is a `phyDat` object it is first converted via
[`PhyDatToMatrix()`](https://ms609.github.io/TreeTools/reference/MatrixToPhyDat.md)
with `ambigNA = TRUE, inappNA = TRUE`, so that fully-ambiguous and
inapplicable rows become `NA_integer_` and only true partial
polymorphisms are subject to the `polymorphism` rule.

## See also

Other phylogenetic matrix conversion functions:
[`Decompose()`](https://ms609.github.io/TreeTools/reference/Decompose.md),
[`MatrixToPhyDat()`](https://ms609.github.io/TreeTools/reference/MatrixToPhyDat.md),
[`Reweight()`](https://ms609.github.io/TreeTools/reference/Reweight.md),
[`StringToPhyDat()`](https://ms609.github.io/TreeTools/reference/PhyToString.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tokens <- matrix(c("0", "(12)", "1", "?", "-"),
                 nrow = 1,
                 dimnames = list("Taxon_A", paste0("C", 1:5)))
NexusTokensToInteger(tokens)
#>         C1 C2 C3 C4 C5
#> Taxon_A  0 NA  1 NA NA
NexusTokensToInteger(tokens, polymorphism = "first")
#>         C1 C2 C3 C4 C5
#> Taxon_A  0  1  1 NA NA
```
