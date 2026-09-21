# Construct consensus trees

`Consensus()` calculates the majority-rule or strict consensus of a set
of trees, using the cluster-table approach of (Day 1985) .

## Usage

``` r
Consensus(trees, p = 1, check.labels = TRUE, hash = TRUE)
```

## Arguments

- trees:

  List of trees, optionally of class `multiPhylo`.

- p:

  A number from 0.5 to 1 giving the proportion of trees that must
  contain a split for it to be reported in the consensus: from `p = 0.5`
  (more than half the trees; the majority-rule consensus) to `p = 1`
  (every tree; the strict consensus, the default).

- check.labels:

  Logical specifying whether to check that all trees have identical
  labels. Defaults to `TRUE`, which is slower.

- hash:

  Logical; if `TRUE` (default), majority/threshold consensus counts
  splits using 128-bit hashing, which is exact with overwhelming
  probability (a collision conflating two distinct splits is vanishingly
  unlikely). Set `hash = FALSE` for a slower but guaranteed-exact count.
  Ignored when `p = 1`, which is always exact.

## Value

`Consensus()` returns an object of class `phylo`, rooted as in the first
entry of `trees`.

## Details

The strict consensus (`p = 1`) compares the clusters of the first tree
against every other tree in linear time. The majority-rule and threshold
consensus (`0.5 <= p < 1`) instead count the frequency of every split
across all trees in a single pass and retain those occurring in a
proportion `p` or more of trees (i.e. in at least
`ceiling(p * length(trees))` trees); this runs in time linear in the
number of trees, after (Jansson et al. 2016) . The majority threshold
`p = 0.5` is strict: a split is retained only if it occurs in *more*
than half the trees, so that two conflicting splits can never both be
reported. By default the count uses a 128-bit hash, whose results are
exact with overwhelming probability; set `hash = FALSE` for a slower but
guaranteed-exact count.

## References

Day WHE (1985). “Optimal algorithms for comparing trees with labeled
leaves.” *Journal of Classification*, **2**(1), 7–28.
[doi:10.1007/BF01908061](https://doi.org/10.1007/BF01908061) .  
  
Jansson J, Shen C, Sung W (2016). “Improved algorithms for constructing
consensus trees.” *Journal of the ACM*, **63**(3), 28:1–28:24.
[doi:10.1145/2898436](https://doi.org/10.1145/2898436) .

## See also

- [ConsTree](https://constree.github.io/) implements other consensus
  tree algorithms.

- [Rogue](https://ms609.github.io/Rogue/) increases the resolution of
  consensus trees by dropping wildcard taxa.

- [`TreeDist::ConsensusInfo()`](https://ms609.github.io/TreeDist/reference/TreeInfo.html)
  calculates the information content of a consensus tree.

Other consensus tree functions:
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/dev/reference/ConsensusWithout.md),
[`RoguePlot()`](https://ms609.github.io/TreeTools/dev/reference/RoguePlot.md)

Other tree characterization functions:
[`CladisticInfo()`](https://ms609.github.io/TreeTools/dev/reference/CladisticInfo.md),
[`J1Index()`](https://ms609.github.io/TreeTools/dev/reference/J1Index.md),
[`Stemwardness`](https://ms609.github.io/TreeTools/dev/reference/Stemwardness.md),
[`TotalCopheneticIndex()`](https://ms609.github.io/TreeTools/dev/reference/TotalCopheneticIndex.md)

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
