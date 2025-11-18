# Loading phylogenetic data into R

It can be a bit fiddly to get a phylogenetic dataset into R,
particularly if you are not used to working with files in the Nexus
format.

First off, make sure that you are comfortable [telling R where to find a
file](https://ms609.github.io/TreeTools/dev/articles/filesystem-navigation.md).

Then you are ready to load raw data:

## Load raw data

### From an Excel spreadsheet

If your data is in an Excel spreadsheet, one way to load it into R is
using the ‘[readxl](https://readxl.tidyverse.org/)’ package. First
you’ll have to install it:

``` r
install.packages("readxl") # You only need to do this once
```

Then you should prepare your Excel spreadsheet such that each row
corresponds to a taxon, and each column to a character.

Then you can read the data from the Excel file by telling R which sheet,
rows and columns contain your data:

``` r
library("readxl")
raw_data <- as.matrix(read_excel(
  filename,
  sheet = 1,           # Loads sheet number 1 from the excel file
  range = "B1:AA21",   # Extracts columns B to AA, rows 1 to 21
  # Note that the first row is interpreted as column (character) names
  col_types = "text"   # Read all columns as character strings
))

# Read row (taxon) names from column A
# Again, the first cell will be interpreted as a column name
taxon_names <- unlist(read_excel(filename, sheet = 1, range = "A1:A21"))

rownames(raw_data) <- taxon_names
```

### From a text or CSV (comma separated values) file

Characters can be read from a text file in a similar manner to Excel.
You may need to adjust the R commands to match the particular format of
your input file.

``` r
raw_data <- read.table(
  filename,            # Path to your input file
  sep = ",",           # What character separates columns?
  header = TRUE,       # Does the data contain a header row?
  row.names = 1,       # Which column contains the row names?
  na.strings = "", 
  stringsAsFactors = FALSE
)
```

### From a Nexus file

TreeTools contains an inbuilt Nexus parser:

``` r
raw_data <- ReadCharacters(filename)
# Or, to go straight to PhyDat format:
as_phydat <- ReadAsPhyDat(filename)
```

This will extract character names and codings from a dataset. It’s been
written to work with datasets downloaded from
[MorphoBank](https://morphobank.org/), but my aim is for this function
to handle most valid (and many invalid) NEXUS files. If you find a file
that this function can’t handle, please [let me
know](https://github.com/ms609/TreeTools/issues/new?title=Unsupported+NEXUS+file&body=Please%20attach%20the%20problematic%20NEXUS%20file%20and%20describe%20the%20issue&labels=enhancement)
and I’ll try to fix it.

In the meantime, alternative Nexus parsers are available: try

``` r
raw_data <- ape::read.nexus.data(filename)
```

Non-standard elements of a Nexus file might be beyond the capabilities
of ape’s parser. In particular, you will need to replace spaces in taxon
names with an underscore, and to arrange all data into a single block
starting `BEGIN DATA`. You’ll need to strip out comments, character
definitions and separate taxon blocks.

The function `readNexus` in package `phylobase` uses the NCL library and
promises to be more powerful, but I’ve not been able to get it to work.

### From a TNT file

A TNT format dataset downloaded from
[MorphoBank](https://morphobank.org/) can be parsed with
`ReadTntCharacters`, which might also handle other TNT-compatible files.
If there’s a file that’s not being read correctly, please [let me
know](https://github.com/ms609/TreeTools/issues/new?title=Unsupported+TNT+file&body=Please%20attach%20the%20problematic%20TNT%20file%20and%20describe%20the%20issue&labels=enhancement)
and I’ll try to fix it.

``` r
raw_data <- ReadTntCharacters(filename)
# Or, to go straight to PhyDat format:
my_data <- ReadTntAsPhyDat(filename)
```

## Processing raw data

Next, we need the raw data in the R-friendly `phyDat` format. If you’ve
used the `ReadAsPhyDat` or `ReadTntAsPhyDat` functions, then you can
skip this step – you’re already there.

Otherwise, you can try

``` r
my_data <- PhyDat(raw_data)
```

or if that doesn’t work,

``` r
my_data <- MatrixToPhyDat(raw_data)
```

These functions are pretty robust, but might return an error when they
encounter an unexpected dataset format – if they don’t work on your
dataset, please  
[let me
know](https://github.com/ms609/TreeTools/issues/new?title=data+to+PhyDat+conversion+issue&body=Please%20describe+the+problem+here:&labels=bug).

Failing that, you can enlist the help of the ‘phangorn’ package:

``` r
install.packages("phangorn")
library("phangorn")
my_data <- phyDat(raw_data, type = "USER", levels = c(0:9, "-"))
```

`type="USER"` tells the parser to expect morphological data.

The `levels` parameter simply lists all the states that any character
might take. `0:9` includes all the integer digits from 0 to 9. If you
have inapplicable data in your matrix, you should list `-` as a separate
level as it represents an additional state (as handled by the Morphy
implementation of (Brazeau, Guillerme, & Smith, 2019)). If you have more
complicated ambiguities, you may need to use a contrast matrix to decode
your matrix.

A contrast matrix translates the tokens used in your dataset to the
character states to which they correspond: for example decoding ‘A’ to
{01}. For more details, see the ‘phangorn-specials’ vignette in the
phangorn package, accessible by typing ‘?phangorn’ in the R prompt and
navigating to index \> package vignettes.

``` r
contrast.matrix <- matrix(data = c(
# 0 1 -  # Each column corresponds to a character-state
  1, 0, 0, # Each row corresponds to a token, here 0, denoting the 
           # character-state set {0} 
  0, 1, 0, # 1 | {1}
  0, 0, 1, # - | {-}
  1, 1, 0, # A | {01}
  1, 1, 0, # + | {01}
  1, 1, 1  # ? | {01-}
), ncol = 3, # ncol corresponds to the number of columns in the matrix
byrow = TRUE)
dimnames(contrast.matrix) <- list(
  c(0, 1, "-", "A", "+", "?"), # A list of the tokens corresponding to each row
                               # in the contrast matrix
  c(0, 1, "-") # A list of the character-states corresponding to the columns 
               # in the contrast matrix
)

contrast.matrix
```

    ##   0 1 -
    ## 0 1 0 0
    ## 1 0 1 0
    ## - 0 0 1
    ## A 1 1 0
    ## + 1 1 0
    ## ? 1 1 1

If you need to use a contrast matrix, convert the data using

``` r
my.phyDat <- phyDat(my.data, type = "USER", contrast = contrast.matrix)
```

## Store processed data

To simplify subsequent analysis, or to allow the data to be read into
other contexts (e.g. via the
[`TreeSearch::EasyTrees()`](https://ms609.github.io/TreeSearch/articles/tree-search.html)
interface), you could save your processed data in Nexus format:

``` r
ape::write.nexus.data(my.phyDat, "my_data.nex", format = "standard")
```

## What next?

You might want to:

- [Load a phylogenetic
  tree](https://ms609.github.io/TreeTools/dev/articles/load-trees.md)
  into R.

- Conduct parsimony search using Brazeau, Guillerme & Smith’s (2019)
  [approach to inapplicable
  data](https://ms609.github.io/TreeSearch/articles/tree-search.html),
  or using [Profile
  parsimony](https://ms609.github.io/TreeSearch/articles/profile.html)
  (Faith & Trueman, 2001).

## References

Brazeau, M. D., Guillerme, T., & Smith, M. R. (2019). An algorithm for
morphological phylogenetic analysis with inapplicable data. *Systematic
Biology*, *68*, 619–631. doi:
[10.1093/sysbio/syy083](https://doi.org/10.1093/sysbio/syy083)

Faith, D. P., & Trueman, J. W. H. (2001). Towards an inclusive
philosophy for phylogenetic inference. *Systematic Biology*, *50*(3),
331–350. doi:
[10.1080/10635150118627](https://doi.org/10.1080/10635150118627)
