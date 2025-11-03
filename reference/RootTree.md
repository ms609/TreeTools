# Root or unroot a phylogenetic tree

`RootTree()` roots a tree on the smallest clade containing the specified
tips; `RootOnNode()` roots a tree on a specified internal node;
`UnrootTree()` collapses a root node, without the undefined behaviour
encountered when using
[`ape::unroot()`](https://rdrr.io/pkg/ape/man/root.html) on trees in
preorder.

## Usage

``` r
RootTree(tree, outgroupTips, fallback = NULL)

RootOnNode(tree, node, resolveRoot = FALSE)

UnrootTree(tree)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html),
  or a list of trees of class `list` or `multiPhylo`.

- outgroupTips:

  Vector of type character, integer or logical, specifying the names or
  indices of the tips to include in the outgroup. If `outgroupTips` is a
  of type character, and a tree contains multiple tips with a matching
  label, the first will be used.

- fallback:

  Vector corresponding to `outgroupTips` determining behaviour when
  `outgroupTips` do not root a tree.

  Where the smallest clade that contains `outgroupTips` includes all
  taxa, `RootTree()` will not change the topology of a tree. If
  `fallback = NULL`, `RootTree()` will return `tree`. Otherwise, taxa
  will be excluded from `outgroupTips` in the sequence specified
  (`fallback[1]` first) by until the clade containing all outgroup tips
  is smaller than `tree`.

- node:

  Integer specifying node (internal or tip) to set as the root.

- resolveRoot:

  Logical specifying whether to resolve the root node.

## Value

`RootTree()` returns a tree of class `phylo`, rooted on the smallest
clade that contains the specified tips, with edges and nodes numbered in
preorder. Node labels are not retained.

`RootOnNode()` returns a tree of class `phylo`, rooted on the requested
`node` and ordered in
[`Preorder`](https://ms609.github.io/TreeTools/reference/Reorder.md).

`UnrootTree()` returns `tree`, in preorder, having collapsed the first
child of the root node in each tree.

## See also

- [`ape::root()`](https://rdrr.io/pkg/ape/man/root.html)

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
[`SortTree()`](https://ms609.github.io/TreeTools/reference/SortTree.md),
[`Subtree()`](https://ms609.github.io/TreeTools/reference/Subtree.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/reference/TipTimedTree.md),
[`TrivialTree`](https://ms609.github.io/TreeTools/reference/TrivialTree.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tree <- PectinateTree(8)
plot(tree)
ape::nodelabels()


plot(RootTree(tree, c("t6", "t7")))


plot(RootOnNode(tree, 12))

plot(RootOnNode(tree, 2))

```
