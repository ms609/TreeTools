# Relative length of internal branches

Treeness (also termed stemminess) is the proportion of a tree's length
found on internal branches (Lanyon 1988) . Insofar as external branches
do not contain phylogenetic (grouping) signal, trees with a high
treeness can be interpreted as containing a higher signal:noise ratio
(Phillips and Penny 2003-08) .

## Usage

``` r
Treeness(tree)

Stemminess(tree)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html),
  or a list of trees of class `list` or `multiPhylo`.

## Value

`Treeness()` returns a numeric vector reporting the treeness of each
`tree`.

## References

Lanyon SM (1988). “The Stochastic Mode of Molecular Evolution: What
Consequences for Systematic Investigations?” *The Auk*, **105**(3),
565–573.
[doi:10.1093/auk/105.3.565](https://doi.org/10.1093/auk/105.3.565) .  
  
Phillips MJ, Penny D (2003-08). “The Root of the Mammalian Tree Inferred
from Whole Mitochondrial Genomes.” *Molecular Phylogenetics and
Evolution*, **28**(2), 171–185.
[doi:10.1016/S1055-7903(03)00057-5](https://doi.org/10.1016/S1055-7903%2803%2900057-5)
.

## See also

Other tree properties:
[`Cherries()`](https://ms609.github.io/TreeTools/dev/reference/Cherries.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/dev/reference/ConsensusWithout.md),
[`LongBranch()`](https://ms609.github.io/TreeTools/dev/reference/LongBranch.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/dev/reference/MatchEdges.md),
[`NSplits()`](https://ms609.github.io/TreeTools/dev/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/dev/reference/NTip.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/dev/reference/NodeNumbers.md),
[`PathLengths()`](https://ms609.github.io/TreeTools/dev/reference/PathLengths.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/dev/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/dev/reference/TipLabels.md),
[`TreeIsRooted()`](https://ms609.github.io/TreeTools/dev/reference/TreeIsRooted.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
lowTree <- BalancedTree(6, lengths = c(1, 1, 4, 4, 4, 1, 1, 4, 4, 4))
plot(lowTree)

Treeness(lowTree)
#> [1] 0.1428571
highTree <- BalancedTree(6, lengths = c(6, 6, 1, 1, 1, 6, 6, 1, 1, 1))
plot(highTree)

Treeness(c(lowTree, highTree))
#> [1] 0.1428571 0.8000000
```
