# TreeTools 0.1.4.9002 (development)

- `NodeDepth()`: discriminate shallow from deep nodes.
- Experiments with Clustering Information.

# TreeTools 0.1.4.9001 (development)

- Improve support for unrooted trees in `as.Splits()`

# TreeTools 0.1.4.9000 (development)

- `RootNode()`: Which node is a tree's root?
- `PhylogeneticInfo()`: Calculate the information content of a tree.
- `NodeOrder()`, `NDescendants()`: Calculate edges incident to each node.
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
