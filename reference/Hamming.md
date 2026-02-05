# Hamming distance between taxa in a phylogenetic dataset

The Hamming distance between a pair of taxa is the number of characters
with a different coding, i.e. the smallest number of evolutionary steps
that must have occurred since their common ancestor.

## Usage

``` r
Hamming(
  dataset,
  ratio = TRUE,
  ambig = c("median", "mean", "zero", "one", "na", "nan")
)
```

## Arguments

- dataset:

  Object of class `phyDat`.

- ratio:

  Logical specifying whether to weight distance against maximum
  possible, given that a token that is ambiguous in either of two taxa
  cannot contribute to the total distance between the pair.

- ambig:

  Character specifying value to return when a pair of taxa have a zero
  maximum distance (perhaps due to a preponderance of ambiguous tokens).
  "median", the default, take the median of all other distance values;
  "mean", the mean; "zero" sets to zero; "one" to one; "NA" to
  `NA_integer_`; and "NaN" to `NaN`.

## Value

`Hamming()` returns an object of class `dist` listing the Hamming
distance between each pair of taxa.

## Details

Tokens that contain the inapplicable state are treated as requiring no
steps to transform into any applicable token.

## See also

Used to construct neighbour joining trees in
[`NJTree()`](https://ms609.github.io/TreeTools/reference/NJTree.md).

`dist.hamming()` in the phangorn package provides an alternative
implementation.

Other utility functions:
[`ClusterTable`](https://ms609.github.io/TreeTools/reference/ClusterTable.md),
[`ClusterTable-methods`](https://ms609.github.io/TreeTools/reference/ClusterTable-methods.md),
[`MSTEdges()`](https://ms609.github.io/TreeTools/reference/MSTEdges.md),
[`SampleOne()`](https://ms609.github.io/TreeTools/reference/SampleOne.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/reference/TipTimedTree.md),
[`UnshiftTree()`](https://ms609.github.io/TreeTools/reference/UnshiftTree.md),
[`as.multiPhylo()`](https://ms609.github.io/TreeTools/reference/as.multiPhylo.md),
[`match,phylo,phylo-method`](https://ms609.github.io/TreeTools/reference/match.multiPhylo.md),
[`sapply64()`](https://ms609.github.io/TreeTools/reference/sapply64.md),
[`sort.multiPhylo()`](https://ms609.github.io/TreeTools/reference/sort.multiPhylo.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tokens <- matrix(c(0, 0, "0", 0, "?",
                   0, 0, "1", 0, 1,
                   0, 0, "1", 0, 1,
                   0, 0, "2", 0, 1,
                   1, 1, "-", "?", 0,
                   1, 1, "2", 1, "{01}"),
                   nrow = 6, ncol = 5, byrow = TRUE,
                   dimnames = list(
                     paste0("Taxon_", LETTERS[1:6]),
                     paste0("Char_", 1:5)))

dataset <- MatrixToPhyDat(tokens)
Hamming(dataset)
#>         Taxon_A Taxon_B Taxon_C Taxon_D Taxon_E
#> Taxon_B    0.25                                
#> Taxon_C    0.25    0.00                        
#> Taxon_D    0.25    0.20    0.20                
#> Taxon_E    1.00    1.00    1.00    1.00        
#> Taxon_F    1.00    0.80    0.80    0.60    0.00
```
