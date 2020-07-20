# TreeTools 1.1.0.9000

 - Documentation improvements.
 - Improve TNT dataset parsing.

# TreeTools 1.1.0

## New functions

 - `RandomTree()`: Draw tree from uniform distribution, instead of via
 `ape::rtree()`.
 - `MakeTreeBinary()`: Uniform equivalent of `ape::multi2di()`.
 - `match.list()` method for lists of splits.
 - `SplitsInBinaryTree()`: How many splits occur in an _n_-leaf binary tree?
 - `vapply64()`, `sapply64()`, `replicate64()`: helper functions when a function
 returns a 64-bit integer.

## Enhancements

 - Use methods for `UnrootTree()`, `RootTree()`, `RootOnNode()` to support
 lists of trees.


# TreeTools 1.0.0

## New functions

- `CladisticInfo()`: Calculate the information content of a tree.
- `RootNode()`: Which node is a tree's root?
- `UnrootTree()`: Safely remove a root node.
- `NodeDepth()`: Discriminate shallow from deep nodes.
- `NodeOrder()`, `NDescendants()`: Count edges incident to each node.
- `CladeSizes()`: Count leaves / nodes descended from each node.
- `ListAncestors()`: List ancestors of a node.
- `LabelSplits()`: Label splits on plotted tree.
- `DropTip()`: Remove tip, handling weird node orders.
- `LeafLabelInterchange()`: Exchange position of _n_ tips.
- `StarTree()`: Generate unresolved tree.
- `TotalCopheneticIndex()` integrated from 'tci' package.

## Deprecations

- `PostorderEdges()`: use `Postorder()` instead.
- `NonDuplicateRoot()`: unused internal function.
- `match.Splits()`: use `match()` instead.
- `in.Splits()`: use `%in%.Splits()` instead.

## Enhancements

- Improve support for unrooted trees in `as.Splits()`.
- Use methods so `Reorder` functions can handle `multiPhylo` objects and edges.
- Handle funny node orders.
- Support continuous characters in `ReadCharacters()`.
- Improve performance of `as.logical.Splits()` and related functions.
- Fail nicely when trees are too large for memory.
- Fix memory leak in `as.Splits()`.
- Various under-the-hood improvements to functions.
- Documentation improvements.


# TreeTools 0.1.4

- Catch hang-inducing bugs in `RootOnNode()`.
- Update `doubleFactorials` cache to fix `as.integer()` rounding error.
- Support unrooted trees in `AddTipEverywhere()`.
- Documentation improvements.


# TreeTools 0.1.3

- `RootOnNode()`: Quickly root a tree on a specified node.
- Improve portability of C++ code.


# TreeTools 0.1.2
 
- `as.Newick`: Fast conversion to Newick format.
- `as.TreeNumber`: Tree shape enumeration.


# TreeTools 0.1.1
 
- Add functions to translate trees to mixed base integers.
- `RenumberTips` can extract tip order from `phylo` and `Splits` objects.
- Documentation changes to attempt to satisfy CRAN submission requirements.


# TreeTools 0.1.0

- Pre-release version spun out of ['TreeSearch'](https://ms609.github.io/TreeSearch)
  package.  Some functionality is subject to change.
