# Reduced consensus, omitting specified taxa

`ConsensusWithout()` displays a consensus plot with specified taxa
excluded, which can be a useful way to increase the resolution of a
consensus tree when a few wildcard taxa obscure a consistent set of
relationships. `MarkMissing()` adds missing taxa as loose leaves on the
plot.

## Usage

``` r
ConsensusWithout(trees, tip = character(0), ...)

# S3 method for class 'phylo'
ConsensusWithout(trees, tip = character(0), ...)

# S3 method for class 'multiPhylo'
ConsensusWithout(trees, tip = character(0), ...)

# S3 method for class 'list'
ConsensusWithout(trees, tip = character(0), ...)

MarkMissing(tip, position = "bottomleft", ...)
```

## Arguments

- trees:

  A list of phylogenetic trees, of class `multiPhylo` or `list`.

- tip:

  A character vector specifying the names (or numbers) of tips to drop
  (using
  [`ape::drop.tip()`](https://rdrr.io/pkg/ape/man/drop.tip.html)).

- ...:

  Additional parameters to pass on to
  [`ape::consensus()`](https://rdrr.io/pkg/ape/man/consensus.html) or
  [`legend()`](https://rdrr.io/r/graphics/legend.html).

- position:

  Where to plot the missing taxa. See
  [`legend()`](https://rdrr.io/r/graphics/legend.html) for options.

## Value

`ConsensusWithout()` returns a consensus tree (of class `phylo`) without
the excluded taxa.

`MarkMissing()` provides a null return, after plotting the specified
`tip`s as a legend.

## See also

Other tree manipulation:
[`AddTip()`](https://ms609.github.io/TreeTools/reference/AddTip.md),
[`CollapseNode()`](https://ms609.github.io/TreeTools/reference/CollapseNode.md),
[`DropTip()`](https://ms609.github.io/TreeTools/reference/DropTip.md),
[`ImposeConstraint()`](https://ms609.github.io/TreeTools/reference/ImposeConstraint.md),
[`KeptPaths()`](https://ms609.github.io/TreeTools/reference/KeptPaths.md),
[`KeptVerts()`](https://ms609.github.io/TreeTools/reference/KeptVerts.md),
[`LeafLabelInterchange()`](https://ms609.github.io/TreeTools/reference/LeafLabelInterchange.md),
[`MakeTreeBinary()`](https://ms609.github.io/TreeTools/reference/MakeTreeBinary.md),
[`Renumber()`](https://ms609.github.io/TreeTools/reference/Renumber.md),
[`RenumberTips()`](https://ms609.github.io/TreeTools/reference/RenumberTips.md),
[`RenumberTree()`](https://ms609.github.io/TreeTools/reference/Reorder.md),
[`RootTree()`](https://ms609.github.io/TreeTools/reference/RootTree.md),
[`SortTree()`](https://ms609.github.io/TreeTools/reference/SortTree.md),
[`Subtree()`](https://ms609.github.io/TreeTools/reference/Subtree.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/reference/TipTimedTree.md),
[`TrivialTree`](https://ms609.github.io/TreeTools/reference/TrivialTree.md)

Other tree properties:
[`Cherries()`](https://ms609.github.io/TreeTools/reference/Cherries.md),
[`LongBranch()`](https://ms609.github.io/TreeTools/reference/LongBranch.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/reference/MatchEdges.md),
[`NSplits()`](https://ms609.github.io/TreeTools/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/reference/NTip.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/reference/NodeNumbers.md),
[`PathLengths()`](https://ms609.github.io/TreeTools/reference/PathLengths.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/reference/TipLabels.md),
[`TreeIsRooted()`](https://ms609.github.io/TreeTools/reference/TreeIsRooted.md),
[`Treeness()`](https://ms609.github.io/TreeTools/reference/Treeness.md)

Other consensus tree functions:
[`Consensus()`](https://ms609.github.io/TreeTools/reference/Consensus.md),
[`RoguePlot()`](https://ms609.github.io/TreeTools/reference/RoguePlot.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
oldPar <- par(mfrow = c(1, 2), mar = rep(0.5, 4))

# Two trees differing only in placement of tip 2:
trees <- as.phylo(c(0, 53), 6)
plot(trees[[1]])
plot(trees[[2]])


# Strict consensus (left panel) lacks resolution:
plot(ape::consensus(trees))

# But omitting tip two (right panel) reveals shared structure in common:
plot(ConsensusWithout(trees, "t2"))
MarkMissing("t2")


par(oldPar)
```
