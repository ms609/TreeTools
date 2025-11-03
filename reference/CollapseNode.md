# Collapse nodes on a phylogenetic tree

Collapses specified nodes or edges on a phylogenetic tree, resulting in
polytomies.

## Usage

``` r
CollapseNode(tree, nodes)

# S3 method for class 'phylo'
CollapseNode(tree, nodes)

CollapseEdge(tree, edges)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- nodes, edges:

  Integer vector specifying the nodes or edges in the tree to be
  dropped. (Use
  [`nodelabels()`](https://rdrr.io/pkg/ape/man/nodelabels.html) or
  [`edgelabels()`](https://rdrr.io/pkg/ape/man/nodelabels.html) to view
  numbers on a plotted tree.)

## Value

`CollapseNode()` and `CollapseEdge()` return a tree of class `phylo`,
corresponding to `tree` with the specified nodes or edges collapsed. The
length of each dropped edge will (naively) be added to each descendant
edge.

## See also

Other tree manipulation:
[`AddTip()`](https://ms609.github.io/TreeTools/reference/AddTip.md),
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
[`TipTimedTree()`](https://ms609.github.io/TreeTools/reference/TipTimedTree.md),
[`TrivialTree`](https://ms609.github.io/TreeTools/reference/TrivialTree.md)

## Author

Martin R. Smith

## Examples

``` r
oldPar <- par(mfrow = c(3, 1), mar = rep(0.5, 4))

tree <- as.phylo(898, 7)
tree$edge.length <- 11:22
plot(tree)
nodelabels()
edgelabels()
edgelabels(round(tree$edge.length, 2),
           cex = 0.6, frame = "n", adj = c(1, -1))

# Collapse by node number
newTree <- CollapseNode(tree, c(12, 13))
plot(newTree)
nodelabels()
edgelabels(round(newTree$edge.length, 2),
           cex = 0.6, frame = "n", adj = c(1, -1))

# Collapse by edge number
newTree <- CollapseEdge(tree, c(2, 4))
plot(newTree)


par(oldPar)
```
