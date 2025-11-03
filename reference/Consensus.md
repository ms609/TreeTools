# Construct consensus trees

`Consensus()` calculates the consensus of a set of trees, using the
algorithm of (Day 1985) .

## Usage

``` r
Consensus(trees, p = 1, check.labels = TRUE)
```

## Arguments

- trees:

  List of trees, optionally of class `multiPhylo`.

- p:

  Proportion of trees that must contain a split for it to be reported in
  the consensus. `p = 0.5` gives the majority-rule consensus; `p = 1`
  (the default) gives the strict consensus.

- check.labels:

  Logical specifying whether to check that all trees have identical
  labels. Defaults to `TRUE`, which is slower.

## Value

`Consensus()` returns an object of class `phylo`, rooted as in the first
entry of `trees`.

## References

Day WHE (1985). “Optimal algorithms for comparing trees with labeled
leaves.” *Journal of Classification*, **2**(1), 7–28.
[doi:10.1007/BF01908061](https://doi.org/10.1007/BF01908061) .

## See also

[`TreeDist::ConsensusInfo()`](https://ms609.github.io/TreeDist/reference/TreeInfo.html)
calculates the information content of a consensus tree.

Other consensus tree functions:
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/reference/ConsensusWithout.md),
[`RoguePlot()`](https://ms609.github.io/TreeTools/reference/RoguePlot.md)

Other tree characterization functions:
[`CladisticInfo()`](https://ms609.github.io/TreeTools/reference/CladisticInfo.md),
[`J1Index()`](https://ms609.github.io/TreeTools/reference/J1Index.md),
[`Stemwardness`](https://ms609.github.io/TreeTools/reference/Stemwardness.md),
[`TotalCopheneticIndex()`](https://ms609.github.io/TreeTools/reference/TotalCopheneticIndex.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
Consensus(as.phylo(0:2, 8))
#> 
#> Phylogenetic tree with 8 tips and 6 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.
```
