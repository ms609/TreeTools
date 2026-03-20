# Extract a subtree

`Subtree()` safely extracts a clade from a phylogenetic tree.

## Usage

``` r
Subtree(tree, node)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html),
  with internal numbering in cladewise order (use
  [`Preorder`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md)`(tree)`
  or (slower)
  [`Cladewise`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md)`(tree)`)
  .

- node:

  The number of the node at the base of the clade to be extracted.

## Value

`Subtree()` returns a tree of class `phylo` that represents a clade
extracted from the original tree.

## Details

Modified from the ape function
[`extract.clade`](https://rdrr.io/pkg/ape/man/drop.tip.html), which
sometimes behaves unpredictably. Unlike extract.clade, this function
supports the extraction of "clades" that constitute a single tip.

## See also

Other tree manipulation:
[`AddTip()`](https://ms609.github.io/TreeTools/dev/reference/AddTip.md),
[`CollapseNode()`](https://ms609.github.io/TreeTools/dev/reference/CollapseNode.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/dev/reference/ConsensusWithout.md),
[`DropTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md),
[`ImposeConstraint()`](https://ms609.github.io/TreeTools/dev/reference/ImposeConstraint.md),
[`KeptPaths()`](https://ms609.github.io/TreeTools/dev/reference/KeptPaths.md),
[`KeptVerts()`](https://ms609.github.io/TreeTools/dev/reference/KeptVerts.md),
[`LeafLabelInterchange()`](https://ms609.github.io/TreeTools/dev/reference/LeafLabelInterchange.md),
[`MakeTreeBinary()`](https://ms609.github.io/TreeTools/dev/reference/MakeTreeBinary.md),
[`Renumber()`](https://ms609.github.io/TreeTools/dev/reference/Renumber.md),
[`RenumberTips()`](https://ms609.github.io/TreeTools/dev/reference/RenumberTips.md),
[`RenumberTree()`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md),
[`RootTree()`](https://ms609.github.io/TreeTools/dev/reference/RootTree.md),
[`SortTree()`](https://ms609.github.io/TreeTools/dev/reference/SortTree.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/dev/reference/TipTimedTree.md),
[`TrivialTree`](https://ms609.github.io/TreeTools/dev/reference/TrivialTree.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tree <- Preorder(BalancedTree(8))
plot(tree)
ape::nodelabels()
ape::nodelabels(13, 13, bg="yellow")


plot(Subtree(tree, 13))

```
