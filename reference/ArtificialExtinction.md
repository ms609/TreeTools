# Artificial Extinction

Remove tokens that do not occur in a fossil "template" taxon from a
living taxon, to simulate the process of fossilization in removing data
from a phylogenetic dataset.

## Usage

``` r
ArtificialExtinction(
  dataset,
  subject,
  template,
  replaceAmbiguous = "ambig",
  replaceCoded = "original",
  replaceAll = TRUE,
  sampleFrom = NULL
)

# S3 method for class 'matrix'
ArtificialExtinction(
  dataset,
  subject,
  template,
  replaceAmbiguous = "ambig",
  replaceCoded = "original",
  replaceAll = TRUE,
  sampleFrom = NULL
)

# S3 method for class 'phyDat'
ArtificialExtinction(
  dataset,
  subject,
  template,
  replaceAmbiguous = "ambig",
  replaceCoded = "original",
  replaceAll = TRUE,
  sampleFrom = NULL
)

ArtEx(
  dataset,
  subject,
  template,
  replaceAmbiguous = "ambig",
  replaceCoded = "original",
  replaceAll = TRUE,
  sampleFrom = NULL
)
```

## Arguments

- dataset:

  Phylogenetic dataset of class `phyDat` or `matrix`.

- subject:

  Vector identifying subject taxa, by name or index.

- template:

  Character or integer identifying taxon to use as a template.

- replaceAmbiguous, replaceCoded:

  Character specifying whether tokens that are ambiguous (`?`) or coded
  (not `?`) in the fossil template should be replaced with:

  - `original`: Their original value; i.e. no change;

  - `ambiguous`: The ambiguous token, `?`;

  - `binary`: The tokens `0` or `1`, with equal probability;

  - `uniform`: One of the tokens present in `sampleFrom`, with equal
    probability;

  - `sample`: One of the tokens present in `sampleFrom`, sampled
    according to their frequency.

- replaceAll:

  Logical: if `TRUE`, replace all tokens in a subject; if `FALSE`, leave
  any ambiguous tokens (`?`) ambiguous.

- sampleFrom:

  Vector identifying a subset of characters from which to sample
  replacement tokens. If `NULL`, replacement tokens will be sampled from
  the initial states of all taxa not used as a template (including the
  subjects).

## Value

A dataset with the same class as `dataset` in which entries that are
ambiguous in `template` are made ambiguous in `subject`.

## Details

Further details are provided in Asher and Smith (2022) .

Note: this simple implementation does not account for character
contingency, e.g. characters whose absence imposes inapplicable or
absent tokens on dependent characters.

## References

Asher R, Smith MR (2022). “Phylogenetic signal and bias in
paleontology.” *Systematic Biology*, **71**(4), 986–1008.
[doi:10.1093/sysbio/syab072](https://doi.org/10.1093/sysbio/syab072) .

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
set.seed(1)
dataset <- matrix(c(sample(0:2, 4 * 8, TRUE),
                    "0", "0", rep("?", 6)), nrow = 5,
                    dimnames = list(c(LETTERS[1:4], "FOSSIL"),
                                    paste("char", 1:8)), byrow = TRUE)
artex <- ArtificialExtinction(dataset, c("A", "C"), "FOSSIL")
```
