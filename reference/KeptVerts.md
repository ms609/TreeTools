# Identify vertices retained when leaves are dropped

Identify vertices retained when leaves are dropped

## Usage

``` r
KeptVerts(tree, keptTips, tipLabels = TipLabels(tree))

# S3 method for class 'phylo'
KeptVerts(tree, keptTips, tipLabels = TipLabels(tree))

# S3 method for class 'numeric'
KeptVerts(tree, keptTips, tipLabels = TipLabels(tree))
```

## Arguments

- tree:

  Original tree of class `phylo`, in
  [`Preorder`](https://ms609.github.io/TreeTools/reference/Reorder.md).

- keptTips:

  Either:

  - a logical vector stating whether each leaf should be retained, in a
    sequence corresponding to `tree[["tip.label"]]`; or

  - a character vector listing the leaf labels to retain; or

  - a numeric vector listing the indices of leaves to retain.

- tipLabels:

  Optional character vector naming the leaves of `tree`, if `keptTips`
  is not logical. Inferred from `tree` if unspecified.

## See also

Other tree manipulation:
[`AddTip()`](https://ms609.github.io/TreeTools/reference/AddTip.md),
[`CollapseNode()`](https://ms609.github.io/TreeTools/reference/CollapseNode.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/reference/ConsensusWithout.md),
[`DropTip()`](https://ms609.github.io/TreeTools/reference/DropTip.md),
[`ImposeConstraint()`](https://ms609.github.io/TreeTools/reference/ImposeConstraint.md),
[`KeptPaths()`](https://ms609.github.io/TreeTools/reference/KeptPaths.md),
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

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
master <- BalancedTree(12)
master <- Preorder(master) # Nodes must be listed in Preorder sequence
plot(master)
nodelabels()


allTips <- master[["tip.label"]]
keptTips <- sample(allTips, 8)
plot(KeepTip(master, keptTips))

kept <- KeptVerts(master, allTips %in% keptTips)

map <- which(kept)
# Node `i` in the reduced tree corresponds to node `map[i]` in the original.
```
