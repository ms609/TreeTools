# Extract taxa from a matrix block

Extract leaf labels and character states from a Nexus-formatted matrix.

## Usage

``` r
ExtractTaxa(matrixLines, character_num = NULL, continuous = FALSE)

NexusTokens(tokens, character_num = NULL)
```

## Arguments

- matrixLines:

  Character vector containing lines of a file that include a
  phylogenetic matrix. See
  [`ReadCharacters()`](https://ms609.github.io/TreeTools/dev/reference/ReadCharacters.md)
  for expected format.

- character_num:

  Index of character(s) to return. `NULL`, the default, returns all
  characters.

- continuous:

  Logical specifying whether characters are continuous. Treated as
  discrete if `FALSE`.

- tokens:

  Vector of character strings corresponding to phylogenetic tokens.

## Value

`ExtractTaxa()` returns a matrix with *n* rows, each named for the
relevant taxon, and *c* columns, each corresponding to the respective
character specified in `character_num`.

`NexusTokens()` returns a character vector in which each entry
corresponds to the states of a phylogenetic character, or a list
containing an error message if input is invalid.

## Examples

``` r
fileName <- paste0(system.file(package = "TreeTools"),
                   "/extdata/input/dataset.nex")
matrixLines <- readLines(fileName)[6:11]
ExtractTaxa(matrixLines)
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#> ____1 "'"  "C"  "h"  "a"  "r"  "a"  "c"  "t"  "e"  "r"   "o"   "n"   "e"  
#> ____2 "C"  "h"  "a"  "r"  "a"  "c"  "t"  "e"  "r"  "_"   "t"   "w"   "o"  
#> ____3 "'"  "l"  "o"  "t"  "s"  "-"  "o"  "f"  "-"  "p"   "u"   "n"   "c"  
#> ____4 "C"  "h"  "a"  "r"  "a"  "c"  "t"  "e"  "r"  "_"   "n"   "/"   "a"  
#> ____5 "C"  "h"  "a"  "r"  "a"  "c"  "t"  "e"  "r"  "_"   "5"   "/"   "s"  
#> ____6 "C"  "h"  "a"  "r"  "a"  "c"  "t"  "e"  "r"  "_"   "6"   "/"   NA   
#>       [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25]
#> ____1 "'"   "/"   "a"   "b"   "s"   "e"   "n"   "t"   "p"   "r"   "e"   "s"  
#> ____2 "/"   "a"   "b"   "s"   "e"   "n"   "t"   "p"   "r"   "e"   "s"   "e"  
#> ____3 "t"   "u"   "a"   "t"   "i"   "o"   "n"   ","   "a"   "n"   "d"   "\"" 
#> ____4 "_"   "l"   "o"   "n"   "g"   "_"   "d"   "e"   "s"   "c"   "r"   "i"  
#> ____5 "i"   "m"   "p"   "l"   "e"   "m"   "o"   "r"   "e"   "_"   "c"   "o"  
#> ____6 NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA   
#>       [,26] [,27] [,28] [,29]
#> ____1 "e"   "n"   "t"   ","  
#> ____2 "n"   "t"   ","   NA   
#> ____3 "s"   "o"   "o"   "n"  
#> ____4 "p"   "t"   "i"   "o"  
#> ____5 "m"   "p"   "l"   "e"  
#> ____6 NA    NA    NA    NA   

NexusTokens("01[01]-?")
#>      [,1] [,2] [,3]   [,4] [,5]
#> [1,] "0"  "1"  "[01]" "-"  "?" 
```
