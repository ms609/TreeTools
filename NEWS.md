# TreeTools 0.2.0

## New functions

- `CladisticInfo()`: Calculate the information content of a tree.
- `RootNode()`: Which node is a tree's root?
- `NodeDepth()`: Discriminate shallow from deep nodes.
- `NodeOrder()`, `NDescendants()`: Count edges incident to each node.

## Enhancements

- Improve support for unrooted trees in `as.Splits()`.
- Handle funny node orders in `CollapseNode()`.
- Improve performance of `as.logical.Splits()` and related functions.
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
