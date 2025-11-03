# Display time-calibrated tree using tip information only

`TipTimedTree()` plots a phylogenetic tree against time using an *ad
hoc* approach based on dates associated with the leaves. Nodes are dated
to the youngest possible value, plus an additional "buffer" (specified
with `minEdge`) to ensure that branching order is readable.

## Usage

``` r
TipTimedTree(tree, tipAge, minEdge = 1)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- tipAge:

  Numeric vector specifying the age (in units-of-time ago) associated
  with each tip in `tree$tip.label` in turn. Older ages signify earlier
  tips.

- minEdge:

  Minimum length of edge to allow (in units-of-time)

## Value

`TipTimedTree()` returns a tree with edge lengths set based on the ages
of each tip.

## Details

This experimental function is liable to change its behaviour, or to be
deprecated, in coming releases. Please contact the maintainer if you
find it useful, so that a production-ready version can be prioritized.

## See also

Other utility functions:
[`ClusterTable`](https://ms609.github.io/TreeTools/reference/ClusterTable.md),
[`ClusterTable-methods`](https://ms609.github.io/TreeTools/reference/ClusterTable-methods.md),
[`Hamming()`](https://ms609.github.io/TreeTools/reference/Hamming.md),
[`MSTEdges()`](https://ms609.github.io/TreeTools/reference/MSTEdges.md),
[`SampleOne()`](https://ms609.github.io/TreeTools/reference/SampleOne.md),
[`UnshiftTree()`](https://ms609.github.io/TreeTools/reference/UnshiftTree.md),
[`as.multiPhylo()`](https://ms609.github.io/TreeTools/reference/as.multiPhylo.md),
[`match,phylo,phylo-method`](https://ms609.github.io/TreeTools/reference/match.multiPhylo.md),
[`sapply64()`](https://ms609.github.io/TreeTools/reference/sapply64.md),
[`sort.multiPhylo()`](https://ms609.github.io/TreeTools/reference/sort.multiPhylo.md)

Other tree manipulation:
[`AddTip()`](https://ms609.github.io/TreeTools/reference/AddTip.md),
[`CollapseNode()`](https://ms609.github.io/TreeTools/reference/CollapseNode.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/reference/ConsensusWithout.md),
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
[`TrivialTree`](https://ms609.github.io/TreeTools/reference/TrivialTree.md)

## Examples

``` r
tree <- BalancedTree(6)
plot(TipTimedTree(tree, tipAge = 1:6, minEdge = 2))
```
