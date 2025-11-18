# Convert between matrices and `phyDat` objects

`MatrixToPhyDat()` converts a matrix of tokens to a `phyDat` object;
`PhyDatToMatrix()` converts a `phyDat` object to a matrix of tokens.

## Usage

``` r
MatrixToPhyDat(tokens, tipLabels = rownames(tokens))

PhyDatToMatrix(
  dataset,
  ambigNA = FALSE,
  inappNA = ambigNA,
  parentheses = c("{", "}"),
  sep = ""
)
```

## Arguments

- tokens:

  Matrix of tokens, possibly created with
  [`ReadCharacters()`](https://ms609.github.io/TreeTools/dev/reference/ReadCharacters.md)
  or
  [`ReadTntCharacters()`](https://ms609.github.io/TreeTools/dev/reference/ReadCharacters.md).
  Row names should correspond to leaf labels; column names may
  optionally correspond to character labels.

- tipLabels:

  Optionally, labels for leaves; will override rownames, if specified.

- dataset:

  A dataset of class `phyDat`.

- ambigNA, inappNA:

  Logical specifying whether to denote ambiguous / inapplicable
  characters as `NA` values.

- parentheses:

  Character vector specifying style of parentheses with which to enclose
  ambiguous characters. `c("[", "]")` or `"[]"` will render `[01]`.
  `NULL` will use the token specified in the `phyDat` object; but beware
  that this will be treated as a distinct (non-ambiguous) token if
  re-encoding with `PhyDatToMatrix()`.

- sep:

  Character with which to separate ambiguous tokens, e.g. `','` will
  render `[0,1]`.

## Value

`MatrixToPhyDat()` returns an object of class `phyDat`.

`PhyDatToMatrix()` returns a matrix corresponding to the uncompressed
character states within a `phyDat` object.

## See also

Other phylogenetic matrix conversion functions:
[`Decompose()`](https://ms609.github.io/TreeTools/dev/reference/Decompose.md),
[`Reweight()`](https://ms609.github.io/TreeTools/dev/reference/Reweight.md),
[`StringToPhyDat()`](https://ms609.github.io/TreeTools/dev/reference/PhyToString.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tokens <- matrix(c(0, 0, "0", 0, 0,
                   0, 0, "1", 0, 1,
                   0, 0, "1", 0, 1,
                   0, 0, "2", 0, 1,
                   1, 1, "-", 1, 0,
                   1, 1, "2", 1, "{01}"),
                   nrow = 6, ncol = 5, byrow = TRUE,
                   dimnames = list(
                     paste0("Taxon_", LETTERS[1:6]),
                     paste0("Char_", 1:5)))
                   
MatrixToPhyDat(tokens)
#> $Taxon_A
#> [1] 1 1 1
#> 
#> $Taxon_B
#> [1] 1 2 2
#> 
#> $Taxon_C
#> [1] 1 2 2
#> 
#> $Taxon_D
#> [1] 1 3 2
#> 
#> $Taxon_E
#> [1] 2 4 1
#> 
#> $Taxon_F
#> [1] 2 3 5
#> 
#> attr(,"weight")
#> [1] 3 1 1
#> attr(,"nr")
#> [1] 3
#> attr(,"nc")
#> [1] 4
#> attr(,"index")
#> [1] 1 1 2 1 3
#> attr(,"levels")
#> [1] "-" "0" "1" "2"
#> attr(,"allLevels")
#> [1] "0"    "1"    "2"    "-"    "{01}"
#> attr(,"type")
#> [1] "USER"
#> attr(,"contrast")
#>      - 0 1 2
#> [1,] 0 1 0 0
#> [2,] 0 0 1 0
#> [3,] 0 0 0 1
#> [4,] 1 0 0 0
#> [5,] 0 1 1 0
#> attr(,"class")
#> [1] "phyDat"
data("Lobo", package = "TreeTools")
head(PhyDatToMatrix(Lobo.phy)[, 91:93])
#>                        [,1] [,2] [,3]
#> Tubiluchus_Priapulida  "-"  "0"  "0" 
#> Cricocosmia            "-"  "0"  "0" 
#> Aysheaia               "0"  "0"  "1" 
#> Siberion               "?"  "0"  "?" 
#> Onychodictyon_ferox    "0"  "0"  "1" 
#> Onychodictyon_gracilis "1"  "0"  "1" 
```
