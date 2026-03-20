# Convert between strings and `phyDat` objects

`PhyDatToString()` converts a `phyDat` object as a string;
`StringToPhyDat()` converts a string of character data to a `phyDat`
object.

## Usage

``` r
StringToPhyDat(string, tips, byTaxon = TRUE)

StringToPhydat(string, tips, byTaxon = TRUE)

PhyToString(
  phy,
  parentheses = "{",
  collapse = "",
  ps = "",
  useIndex = TRUE,
  byTaxon = TRUE,
  concatenate = TRUE
)

PhyDatToString(
  phy,
  parentheses = "{",
  collapse = "",
  ps = "",
  useIndex = TRUE,
  byTaxon = TRUE,
  concatenate = TRUE
)

PhydatToString(
  phy,
  parentheses = "{",
  collapse = "",
  ps = "",
  useIndex = TRUE,
  byTaxon = TRUE,
  concatenate = TRUE
)
```

## Arguments

- string:

  String of tokens, optionally containing whitespace, with no
  terminating semi-colon.

- tips:

  (Optional) Character vector corresponding to the names (in order) of
  each taxon in the matrix, or an object such as a tree from which tip
  labels can be extracted.

- byTaxon:

  Logical. If `TRUE`, write one taxon followed by the next. If `FALSE`,
  write one character followed by the next.

- phy:

  An object of class `phyDat`.

- parentheses:

  Character specifying format of parentheses with which to surround
  ambiguous tokens. Choose from: `{` (default), `[`, `(`, `<`.

- collapse:

  Character specifying text, perhaps `,`, with which to separate
  multiple tokens within parentheses.

- ps:

  Character specifying text, perhaps `;`, to append to the end of each
  string.

- useIndex:

  Logical (default: `TRUE`) specifying whether to print duplicate
  characters multiple times, as they appeared in the original matrix.

- concatenate:

  Logical specifying whether to concatenate all characters/taxa into a
  single string, or to return a separate string for each entry.

## Value

`StringToPhyDat()` returns an object of class `phyDat`.

`PhyToString()` returns a character vector listing a text representation
of the phylogenetic character state for each taxon in turn.

## See also

Other phylogenetic matrix conversion functions:
[`Decompose()`](https://ms609.github.io/TreeTools/dev/reference/Decompose.md),
[`MatrixToPhyDat()`](https://ms609.github.io/TreeTools/dev/reference/MatrixToPhyDat.md),
[`Reweight()`](https://ms609.github.io/TreeTools/dev/reference/Reweight.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
StringToPhyDat("-?01231230?-", c("Lion", "Gazelle"), byTaxon = TRUE)
#> $Lion
#> [1] 1 3 5 2 4 6
#> 
#> $Gazelle
#> [1] 2 4 6 5 3 1
#> 
#> attr(,"weight")
#> [1] 1 1 1 1 1 1
#> attr(,"nr")
#> [1] 6
#> attr(,"nc")
#> [1] 5
#> attr(,"index")
#> [1] 1 2 3 4 5 6
#> attr(,"levels")
#> [1] "-" "0" "1" "2" "3"
#> attr(,"levels")attr(,".match.hash")
#> <hash table>
#> attr(,"allLevels")
#> [1] "-" "1" "?" "2" "0" "3"
#> attr(,"type")
#> [1] "USER"
#> attr(,"contrast")
#>      - 0 1 2 3
#> [1,] 1 0 0 0 0
#> [2,] 0 0 1 0 0
#> [3,] 1 1 1 1 1
#> [4,] 0 0 0 1 0
#> [5,] 0 1 0 0 0
#> [6,] 0 0 0 0 1
#> attr(,"class")
#> [1] "phyDat"
# encodes the following matrix:
# Lion     -?0123
# Gazelle  1230?-

fileName <- paste0(system.file(package = "TreeTools"),
                   "/extdata/input/dataset.nex")
phyDat <- ReadAsPhyDat(fileName)
PhyToString(phyDat, concatenate = FALSE)
#> [1] "0000000" "0000000" "1111?00" "111??11" "1111?11"
```
