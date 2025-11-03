# Identify taxa with long branches

The long branch (LB) score (Struck 2014) measures the deviation of the
average pairwise patristic distance of a leaf from all other leaves in a
tree, relative to the average leaf-to-leaf distance.

## Usage

``` r
LongBranch(tree)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html),
  or a list of trees of class `list` or `multiPhylo`.

## Value

`LongBranch()` returns a vector giving the long branch score for each
leaf in `tree`, or a list of such vectors if `tree` is a list. Results
are given as raw deviations, without multiplying by 100 as proposed by
Struck (2014) .

## Details

Struck (2014) proposes the standard deviation of LB scores as a measure
of heterogeneity that can be compared between trees; and the upper
quartile of LB scores as "a representative value for the taxa with the
longest branches".

## See also

Other tree properties:
[`Cherries()`](https://ms609.github.io/TreeTools/reference/Cherries.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/reference/ConsensusWithout.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/reference/MatchEdges.md),
[`NSplits()`](https://ms609.github.io/TreeTools/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/reference/NTip.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/reference/NodeNumbers.md),
[`PathLengths()`](https://ms609.github.io/TreeTools/reference/PathLengths.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/reference/TipLabels.md),
[`TreeIsRooted()`](https://ms609.github.io/TreeTools/reference/TreeIsRooted.md),
[`Treeness()`](https://ms609.github.io/TreeTools/reference/Treeness.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tree <- BalancedTree(8, lengths = c(rep(2, 4), 5:7, rep(2, 4), rep(1, 3)))
lb <- LongBranch(tree)
tree$tip.label <- paste(tree$tip.label, signif(lb, 3), sep = ": ")
plot(tree, tip.col = SupportColour((1 - lb) / 2), font = 2)


# Standard deviation of LB scores allows comparison with other trees
sd(lb)
#> [1] 0.2335139
evenLengths <- BalancedTree(8, lengths = jitter(rep(1, 14)))
sd(LongBranch(evenLengths))
#> [1] 0.001698709

# Upper quartile identifies taxa with longest branches
threshold <- quantile(lb, 0.75)
tree$tip.label[lb > threshold]
#> [1] "t3: 0.333" "t4: 0.403"
```
