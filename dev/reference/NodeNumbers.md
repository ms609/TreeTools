# Numeric index of each node in a tree `NodeNumbers()` returns a sequence corresponding to the nodes in a tree

Numeric index of each node in a tree `NodeNumbers()` returns a sequence
corresponding to the nodes in a tree

## Usage

``` r
NodeNumbers(tree, tips = FALSE)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- tips:

  Logical specifying whether to also include the indices of leaves.

## Value

`NodeNumbers()` returns an integer vector corresponding to the indices
of nodes within a tree.

## See also

Other tree properties:
[`Cherries()`](https://ms609.github.io/TreeTools/dev/reference/Cherries.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/dev/reference/ConsensusWithout.md),
[`LongBranch()`](https://ms609.github.io/TreeTools/dev/reference/LongBranch.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/dev/reference/MatchEdges.md),
[`NSplits()`](https://ms609.github.io/TreeTools/dev/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/dev/reference/NTip.md),
[`PathLengths()`](https://ms609.github.io/TreeTools/dev/reference/PathLengths.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/dev/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/dev/reference/TipLabels.md),
[`TreeIsRooted()`](https://ms609.github.io/TreeTools/dev/reference/TreeIsRooted.md),
[`Treeness()`](https://ms609.github.io/TreeTools/dev/reference/Treeness.md)

Other tree navigation:
[`AncestorEdge()`](https://ms609.github.io/TreeTools/dev/reference/AncestorEdge.md),
[`CladeSizes()`](https://ms609.github.io/TreeTools/dev/reference/CladeSizes.md),
[`DescendantEdges()`](https://ms609.github.io/TreeTools/dev/reference/DescendantEdges.md),
[`EdgeAncestry()`](https://ms609.github.io/TreeTools/dev/reference/EdgeAncestry.md),
[`EdgeDistances()`](https://ms609.github.io/TreeTools/dev/reference/EdgeDistances.md),
[`ListAncestors()`](https://ms609.github.io/TreeTools/dev/reference/ListAncestors.md),
[`MRCA()`](https://ms609.github.io/TreeTools/dev/reference/MRCA.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/dev/reference/MatchEdges.md),
[`NDescendants()`](https://ms609.github.io/TreeTools/dev/reference/NDescendants.md),
[`NodeDepth()`](https://ms609.github.io/TreeTools/dev/reference/NodeDepth.md),
[`NodeOrder()`](https://ms609.github.io/TreeTools/dev/reference/NodeOrder.md),
[`RootNode()`](https://ms609.github.io/TreeTools/dev/reference/RootNode.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)
