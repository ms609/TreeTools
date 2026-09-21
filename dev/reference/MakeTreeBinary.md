# Generate binary tree by collapsing polytomies

`MakeTreeBinary()` resolves, at random, all polytomies in a tree or set
of trees, such that all trees compatible with the input topology are
drawn with equal probability. Edge lengths are not yet supported, so are
removed.

## Usage

``` r
MakeTreeBinary(tree)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

## Value

`MakeTreeBinary()` returns a rooted binary tree of class `phylo`,
corresponding to tree uniformly selected from all those compatible with
the input tree topologies.

## See also

Since ape v5.5, this functionality is available through
[`ape::multi2di()`](https://rdrr.io/pkg/ape/man/multi2di.html); previous
versions of "ape" did not return topologies in equal frequencies.
`MakeTreeBinary()` is often somewhat faster;
[`multi2di()`](https://rdrr.io/pkg/ape/man/multi2di.html) retains edge
lengths.

Other tree manipulation:
[`AddTip()`](https://ms609.github.io/TreeTools/dev/reference/AddTip.md),
[`CollapseNode()`](https://ms609.github.io/TreeTools/dev/reference/CollapseNode.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/dev/reference/ConsensusWithout.md),
[`DropTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md),
[`ImposeConstraint()`](https://ms609.github.io/TreeTools/dev/reference/ImposeConstraint.md),
[`KeptPaths()`](https://ms609.github.io/TreeTools/dev/reference/KeptPaths.md),
[`KeptVerts()`](https://ms609.github.io/TreeTools/dev/reference/KeptVerts.md),
[`LeafLabelInterchange()`](https://ms609.github.io/TreeTools/dev/reference/LeafLabelInterchange.md),
[`Renumber()`](https://ms609.github.io/TreeTools/dev/reference/Renumber.md),
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
MakeTreeBinary(CollapseNode(PectinateTree(7), c(9, 11, 13)))
#> 
#> Phylogenetic tree with 7 tips and 6 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.
UnrootTree(MakeTreeBinary(StarTree(5)))
#> 
#> Phylogenetic tree with 5 tips and 3 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5
#> 
#> Unrooted; no branch length.
```
