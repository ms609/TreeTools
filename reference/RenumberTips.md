# Renumber a tree's tips

`RenumberTips(tree, tipOrder)` sorts the tips of a phylogenetic tree
`tree` such that the indices in `tree[["edge"]][, 2]` correspond to the
order of leaves given in `tipOrder`.

## Usage

``` r
RenumberTips(tree, tipOrder)

# S3 method for class 'phylo'
RenumberTips(tree, tipOrder)

# S3 method for class 'multiPhylo'
RenumberTips(tree, tipOrder)

# S3 method for class 'list'
RenumberTips(tree, tipOrder)

# S3 method for class '`NULL`'
RenumberTips(tree, tipOrder)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- tipOrder:

  A character vector containing the values of `tree[["tip.label"]]` in
  the desired sort order, or an object (perhaps of class `phylo` or
  `Splits`) with tip labels.

## Value

`RenumberTips()` returns `tree`, with the tips' internal representation
numbered to match `tipOrder`.

## See also

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
[`RenumberTree()`](https://ms609.github.io/TreeTools/reference/Reorder.md),
[`RootTree()`](https://ms609.github.io/TreeTools/reference/RootTree.md),
[`SortTree()`](https://ms609.github.io/TreeTools/reference/SortTree.md),
[`Subtree()`](https://ms609.github.io/TreeTools/reference/Subtree.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/reference/TipTimedTree.md),
[`TrivialTree`](https://ms609.github.io/TreeTools/reference/TrivialTree.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
data("Lobo") # Loads the phyDat object Lobo.phy
tree <- RandomTree(Lobo.phy)
tree <- RenumberTips(tree, names(Lobo.phy))
```
