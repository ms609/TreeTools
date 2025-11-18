# Drop leaves from tree

`DropTip()` removes specified leaves from a phylogenetic tree,
collapsing incident branches.

## Usage

``` r
DropTip(tree, tip, preorder = TRUE, check = TRUE)

KeepTip(tree, tip, preorder = TRUE, check = TRUE)

# S3 method for class 'phylo'
DropTip(tree, tip, preorder = TRUE, check = TRUE)

# S3 method for class 'phylo'
KeepTip(tree, tip, preorder = TRUE, check = TRUE)

# S3 method for class 'Splits'
KeepTip(tree, tip, preorder = TRUE, check = TRUE)

# S3 method for class 'Splits'
DropTip(tree, tip, preorder, check = TRUE)

DropTipPhylo(tree, tip, preorder = TRUE, check = TRUE)

# S3 method for class 'multiPhylo'
DropTip(tree, tip, preorder = TRUE, check = TRUE)

# S3 method for class 'multiPhylo'
KeepTip(tree, tip, preorder = TRUE, check = TRUE)

# S3 method for class 'list'
DropTip(tree, tip, preorder = TRUE, check = TRUE)

# S3 method for class 'list'
KeepTip(tree, tip, preorder = TRUE, check = TRUE)

# S3 method for class '`NULL`'
DropTip(tree, tip, preorder = TRUE, check = TRUE)

# S3 method for class '`NULL`'
KeepTip(tree, tip, preorder = TRUE, check = TRUE)

KeepTipPreorder(tree, tip)

KeepTipPostorder(tree, tip)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- tip:

  Character vector specifying labels of leaves in tree to be dropped or
  kept, or integer vector specifying the indices of leaves to be dropped
  or kept. Specifying the index of an internal node will drop all
  descendants of that node.

- preorder:

  Logical specifying whether to
  [Preorder](https://ms609.github.io/TreeTools/dev/reference/Reorder.md)
  `tree` before dropping tips. Specifying `FALSE` saves a little time,
  but will result in undefined behaviour if `tree` is not in preorder.

- check:

  Logical specifying whether to check validity of `tip`. If `FALSE` and
  `tip` contains entries that do not correspond to leaves of the tree,
  undefined behaviour may occur.

## Value

`DropTip()` returns a tree of class `phylo`, with the requested leaves
removed. The edges of the tree will be numbered in preorder, but their
sequence may not conform to the conventions of
[`Preorder()`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md).

`KeepTip()` returns `tree` with all leaves not in `tip` removed, in
preorder.

## Details

This function differs from
[`ape::drop.tip()`](https://rdrr.io/pkg/ape/man/drop.tip.html), which
roots unrooted trees, and which can crash when trees' internal numbering
follows unexpected schema.

## Functions

- `DropTipPhylo()`: Direct call to `DropTip.phylo()`, to avoid overhead
  of querying object's class.

- `KeepTipPreorder()`: Faster version with no checks. Does not retain
  labels or edge weights. Edges must be listed in preorder. May crash if
  improper input is specified.

- `KeepTipPostorder()`: Faster version with no checks. Does not retain
  labels or edge weights. Edges must be listed in postorder. May crash
  if improper input is specified.

## See also

Other tree manipulation:
[`AddTip()`](https://ms609.github.io/TreeTools/dev/reference/AddTip.md),
[`CollapseNode()`](https://ms609.github.io/TreeTools/dev/reference/CollapseNode.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/dev/reference/ConsensusWithout.md),
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
[`Subtree()`](https://ms609.github.io/TreeTools/dev/reference/Subtree.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/dev/reference/TipTimedTree.md),
[`TrivialTree`](https://ms609.github.io/TreeTools/dev/reference/TrivialTree.md)

Other split manipulation functions:
[`SplitConsistent()`](https://ms609.github.io/TreeTools/dev/reference/SplitConsistent.md),
[`Subsplit()`](https://ms609.github.io/TreeTools/dev/reference/Subsplit.md),
[`TrivialSplits()`](https://ms609.github.io/TreeTools/dev/reference/TrivialSplits.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tree <- BalancedTree(9)
plot(tree)

plot(DropTip(tree, c("t5", "t6")))


unrooted <- UnrootTree(tree)
plot(unrooted)

plot(DropTip(unrooted, 4:5))


summary(DropTip(as.Splits(tree), 4:5))
#> 4 bipartition splits dividing 7 tips, t1 .. t9
#>      1234567
#>  12  ***....
#>  13  **.....
#>  16  ...**..
#>  17  .....**
#> 
#>  Tip 1: t1    Tip 2: t2   Tip 3: t3   Tip 4: t6   Tip 5: t7  
#>  Tip 6: t8    Tip 7: t9  
```
