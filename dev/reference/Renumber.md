# Renumber a tree's nodes and tips

`Renumber()` numbers the nodes and tips in a tree to conform with the
`phylo` standards.

## Usage

``` r
Renumber(tree)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

## Value

`Renumber()` returns a tree of class `phylo`, numbered in a
[Cladewise](https://ms609.github.io/TreeTools/dev/reference/Reorder.md)
fashion consistent with the expectations of ape functions.

## Details

The ape class `phylo` is not formally defined, but expects trees'
internal representation to conform to certain principles: for example,
nodes should be numbered sequentially, with values increasing away from
the root.

`Renumber()` attempts to reformat any tree into a representation that
will not cause ape functions to produce unwanted results or to crash R.

## See also

[`Preorder()`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md)
provides a faster and simpler alternative, but also rotates nodes.

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
[`RenumberTips()`](https://ms609.github.io/TreeTools/dev/reference/RenumberTips.md),
[`RenumberTree()`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md),
[`RootTree()`](https://ms609.github.io/TreeTools/dev/reference/RootTree.md),
[`SortTree()`](https://ms609.github.io/TreeTools/dev/reference/SortTree.md),
[`Subtree()`](https://ms609.github.io/TreeTools/dev/reference/Subtree.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/dev/reference/TipTimedTree.md),
[`TrivialTree`](https://ms609.github.io/TreeTools/dev/reference/TrivialTree.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tree <- RandomTree(letters[1:10])
Renumber(tree)
#> 
#> Phylogenetic tree with 10 tips and 8 internal nodes.
#> 
#> Tip labels:
#>   d, g, c, i, f, j, ...
#> 
#> Unrooted; no branch length.
```
