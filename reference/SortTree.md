# Sort tree

`SortTree()` sorts each node into a consistent order, so that node
rotation does not obscure similarities between similar trees.

## Usage

``` r
SortTree(tree, how = "cladesize", order = TipLabels(tree))

# S3 method for class 'phylo'
SortTree(tree, how = "cladesize", order = TipLabels(tree))

# S3 method for class 'list'
SortTree(tree, how = "cladesize", order = TipLabels(tree[[1]]))

# S3 method for class 'multiPhylo'
SortTree(tree, how = "cladesize", order = TipLabels(tree[[1]]))
```

## Arguments

- tree:

  One or more trees of class `phylo`, optionally as a list or a
  `multiPhylo` object.

- how:

  Character vector specifying sort method: `"Cladesize"` rotates each
  node such that the larger clade is first, thus appearing lower when
  plotted; `"TipLabels"` rotates nodes such that labels listed sooner in
  `order` are listed first, and thus plot lower.

- order:

  Character vector listing tip labels in sequence they should appear on
  tree. Clades containing a taxon earlier in this list will be listed
  sooner and thus plot lower on a tree. Taxa not listed in `order` will
  be treated as if they were last in the list.

## Value

`SortTree()` returns tree in the format of `tree`, with each node in
each tree sorted

## Details

At each node, clades will be listed in `tree[["edge"]]` in decreasing
size order.

Clades that contain the same number of leaves are sorted in decreasing
order of minimum leaf number, so (2, 3) will occur before (1, 4).

As trees are plotted from "bottom up", the largest clades will "sink" to
the bottom of a plotted tree.

## See also

[`Preorder()`](https://ms609.github.io/TreeTools/reference/Reorder.md)
also rearranges trees into a consistent shape, based on the index of
leaves.

[`sort.multiPhylo()`](https://ms609.github.io/TreeTools/reference/sort.multiPhylo.md)
sorts a list of trees stored as a `multiPhylo` object.

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
[`Subtree()`](https://ms609.github.io/TreeTools/reference/Subtree.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/reference/TipTimedTree.md),
[`TrivialTree`](https://ms609.github.io/TreeTools/reference/TrivialTree.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
messyTree <- as.phylo(10, 6)
plot(messyTree)


sorted <- SortTree(messyTree)
plot(sorted)
ape::nodelabels()
ape::edgelabels()
ape::tiplabels(adj = c(2, 1/3))


plot(SortTree(messyTree, how = "tip"))
```
