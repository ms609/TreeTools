# Read phylogenetic characters from file

Parse a Nexus (Maddison et al. 1997) or TNT (Goloboff et al. 2008) file,
reading character states and names.

## Usage

``` r
ReadCharacters(filepath, character_num = NULL, encoding = "UTF8")

ReadTntCharacters(
  filepath,
  character_num = NULL,
  type = NULL,
  encoding = "UTF8"
)

ReadTNTCharacters(
  filepath,
  character_num = NULL,
  type = NULL,
  encoding = "UTF8"
)

ReadNotes(filepath, encoding = "UTF8")

ReadAsPhyDat(...)

ReadTntAsPhyDat(...)

ReadTNTAsPhyDat(...)

PhyDat(dataset)
```

## Arguments

- filepath:

  character string specifying location of file, or a
  [connection](https://rdrr.io/r/base/connections.html) to the file.

- character_num:

  Index of character(s) to return. `NULL`, the default, returns all
  characters.

- encoding:

  Character encoding of input file.

- type:

  Character vector specifying categories of data to extract from file.
  Setting `type = c("num", "dna")` will return only characters following
  a `&[num]` or `&[dna]` tag in a TNT input file, listing `num`
  character blocks before `dna` characters. Leave as `NULL` (the
  default) to return all characters in their original sequence.

- ...:

  Parameters to pass to `Read[Tnt]Characters()`.

- dataset:

  list of taxa and characters, in the format produced by
  [`read.nexus.data()`](https://rdrr.io/pkg/ape/man/read.nexus.data.html):
  a list of sequences each made of a single character vector, and named
  with the taxon name.

## Value

`ReadCharacters()` and `ReadTNTCharacters()` return a matrix whose row
names correspond to tip labels, and column names correspond to character
labels, with the attribute `state.labels` listing the state labels for
each character; or a list of length one containing a character string
explaining why the function call was unsuccessful.

`ReadAsPhyDat()` and `ReadTntAsPhyDat()` return a `phyDat` object.

`ReadNotes()` returns a list in which each entry corresponds to a single
character, and itself contains a list of with two elements:

1.  A single character object listing any notes associated with the
    character

2.  A named character vector listing the notes associated with each
    taxon for that character, named with the names of each note-bearing
    taxon.

## Details

Tested with matrices downloaded from
[MorphoBank](https://morphobank.org/) (O’Leary and Kaufman 2011) , but
should also work more widely; please
[report](https://github.com/ms609/TreeTools/issues/new?title=Error+parsing+Nexus+file&body=%3C!--Tell+me+more+and+attach+your+file...--%3E)
incompletely or incorrectly parsed files.

Matrices must contain only continuous or only discrete characters;
maximum one matrix per file. Continuous characters will be read as
strings (i.e. base type "character").

The encoding of an input file will be automatically determined by R.
Errors pertaining to an `invalid multibyte string` or
`string invalid at that locale` indicate that R has failed to detect the
appropriate encoding. Either [re-save the
file](https://support.posit.co/hc/en-us/articles/200532197-Character-Encoding-in-the-RStudio-IDE)
in a supported encoding (`UTF-8` is a good choice) or specify the file
encoding (which you can find by, for example, opening in
[Notepad++](https://notepad-plus-plus.org/downloads/) and identifying
the highlighted option in the "Encoding" menu) following the example
below.

## Functions

- `PhyDat()`: A convenient wrapper for phangorn's `phyDat()`, which
  converts a **list** of morphological characters into a `phyDat`
  object. If your morphological characters are in the form of a
  **matrix**, perhaps because they have been read using
  [`read.table()`](https://rdrr.io/r/utils/read.table.html), try
  [`MatrixToPhyDat()`](https://ms609.github.io/TreeTools/reference/MatrixToPhyDat.md)
  instead.

## References

Goloboff PA, Farris JS, Nixon KC (2008). “TNT, a free program for
phylogenetic analysis.” *Cladistics*, **24**(5), 774–786.  
  
Maddison DR, Swofford DL, Maddison WP (1997). “Nexus: an extensible file
format for systematic information.” *Systematic Biology*, **46**,
590–621.
[doi:10.1093/sysbio/46.4.590](https://doi.org/10.1093/sysbio/46.4.590)
.  
  
O’Leary MA, Kaufman S (2011). “MorphoBank: phylophenomics in the
"cloud".” *Cladistics*, **27**(5), 529–537.

## See also

- Convert between matrices and `phyDat` objects:
  [`MatrixToPhyDat()`](https://ms609.github.io/TreeTools/reference/MatrixToPhyDat.md)

- Write characters to TNT-format file:
  [`WriteTntCharacters()`](https://ms609.github.io/TreeTools/reference/WriteTntCharacters.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
fileName <- paste0(system.file(package = "TreeTools"),
                   "/extdata/input/dataset.nex")
ReadCharacters(fileName)
#>         Character one Character two lots-of-punctuation, and "so on"!
#> taxon_a "0"           "0"           "0"                              
#> taxon_b "0"           "0"           "0"                              
#> taxon_c "1"           "1"           "1"                              
#> taxon_d "1"           "1"           "1"                              
#> taxon_e "1"           "1"           "1"                              
#>         Character n Character 5 Character 6 final character
#> taxon_a "0"         "0"         "0"         "0"            
#> taxon_b "0"         "0"         "0"         "0"            
#> taxon_c "1"         "?"         "0"         "0"            
#> taxon_d "?"         "?"         "1"         "1"            
#> taxon_e "1"         "?"         "1"         "1"            
#> attr(,"state.labels")
#> attr(,"state.labels")[[1]]
#> [1] "absent"  "present"
#> 
#> attr(,"state.labels")[[2]]
#> [1] "absent"  "present"
#> 
#> attr(,"state.labels")[[3]]
#> [1] "here"       "there"      "everywhere"
#> 
#> attr(,"state.labels")[[4]]
#> [1] "a long description" "present"           
#> 
#> attr(,"state.labels")[[5]]
#> [1] "simple"                "more complex"          "with (parentheses)"   
#> [4] "more complex, 6 still"
#> 
#> attr(,"state.labels")[[6]]
#> [1] "this one has"   "multiple lines"
#> 
#> attr(,"state.labels")[[7]]
#> [1] "absent"  "present"
#> 

fileName <- paste0(system.file(package = "TreeTools"),
                   "/extdata/tests/continuous.nex")

continuous <- ReadCharacters(fileName, encoding = "UTF8")

# To convert from strings to numbers:
at <- attributes(continuous)
continuous <- suppressWarnings(as.numeric(continuous))
attributes(continuous) <- at
continuous
#>            [,1]  [,2]  [,3]  [,4]  [,5]  [,6]
#> A_taxon   1.111 1.000 1.330 1.444 1.555 1.666
#> B_alienus 2.111 2.222 2.333    NA 2.550 2.666
#> C_andinus 3.111 3.222 3.333 3.444 3.555 3.666
```
