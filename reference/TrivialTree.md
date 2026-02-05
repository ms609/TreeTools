# Generate trivial trees

`SingleTaxonTree()` creates a phylogenetic "tree" that contains a single
taxon. `ZeroTaxonTree()` creates an empty `phylo` object with zero
leaves or edges.

## Usage

``` r
SingleTaxonTree(label = "t1", lengths = NULL)

ZeroTaxonTree()
```

## Arguments

- label:

  a character vector specifying the label of the tip.

- lengths:

  a numeric vector specifying the edge lengths of the tree.

## Value

`SingleTaxonTree()` returns a `phylo` object containing a single tip
with the specified label.

`ZeroTaxonTree()` returns an empty `phylo` object.

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
[`RenumberTips()`](https://ms609.github.io/TreeTools/reference/RenumberTips.md),
[`RenumberTree()`](https://ms609.github.io/TreeTools/reference/Reorder.md),
[`RootTree()`](https://ms609.github.io/TreeTools/reference/RootTree.md),
[`SortTree()`](https://ms609.github.io/TreeTools/reference/SortTree.md),
[`Subtree()`](https://ms609.github.io/TreeTools/reference/Subtree.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/reference/TipTimedTree.md)

Other tree generation functions:
[`ConstrainedNJ()`](https://ms609.github.io/TreeTools/reference/ConstrainedNJ.md),
[`GenerateTree`](https://ms609.github.io/TreeTools/reference/GenerateTree.md),
[`NJTree()`](https://ms609.github.io/TreeTools/reference/NJTree.md),
[`TreeNumber`](https://ms609.github.io/TreeTools/reference/TreeNumber.md)

## Examples

``` r
SingleTaxonTree("Homo_sapiens")
#> 
#> Phylogenetic tree with 1 tip and 1 internal node.
#> 
#> Tip label:
#>   Homo_sapiens
#> 
#> Rooted; no branch length.
plot(SingleTaxonTree("root") + BalancedTree(4))


ZeroTaxonTree()
#> 
#> Phylogenetic tree with 0 tips and 0 internal nodes.
#> 
#> Tip labels:
#>   
#> 
#> Rooted; no branch length.
```
