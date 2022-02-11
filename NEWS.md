# TreeTools 1.6.0.9008 #

- `as.MixedBase()` uniquely represents binary trees as a mixed-base vector.

- `sort.multiPhylo()` sorts lists of trees according to their mixed base
  representation ([#84](https://github.com/ms609/TreeTools/issues/84)).

- `%in%.Splits()` retains names when comparing small splits
  ([#40](https://github.com/ms609/TreeTools/issues/40)).
  
- Improve handling of zero-split trees.

- `DropTip()` no longer adds a root to unrooted trees.

- Improve speed of `DropTip()`, by an order of magnitude in some cases.

- Support edge lengths in `Preorder()`, `RootTree()`, `UnrootTree()` and
  `Postorder()` ([#49](https://github.com/ms609/TreeTools/issues/49),
  [#89](https://github.com/ms609/TreeTools/issues/89)).

- Fix bug when tree is rooted on a discontinuous outgroup.

- `SortTree()` handles weighted and non-binary trees
  ([#25](https://github.com/ms609/TreeTools/issues/25),
  [#25](https://github.com/ms609/TreeTools/issues/49)),
  and gains option to sort by tip labels.

- Use lighter Rcpp headers.

- Small improvements to computational efficiency.

- `SpectrumLegend()` adds gradients to plot legends.

- Remove deprecated function `PostorderEdges()`
  ([#35](https://github.com/ms609/TreeTools/issues/35)).

- Fix bug handling `phyDat` objects in `ArtificialExtinction()`.

# TreeTools 1.6.0

## New functions

- `RoguePlot()` plots the positions of rogue taxa.


## Improvements

- `DropTip()` gains `check` parameter to allow slightly faster operation where
  input is guaranteed to be valid.

- `RandomTree()` gains `nodes` parameter allow the inclusion of polytomies.

- Infer `tips` parameter if missing in `StringToPhyDat()`.

- Remove dependency on "phangorn" (allowing use on R < 4.1)

- Improve parsing of information from nexus files.

- Export `DropTipPhylo()` as wrapper to `DropTip.phylo()`.


# TreeTools 1.5.1

- `PhyDatToMatrix()` optionally encodes ambiguous / inapplicable tokens as `NA`.

- Implement `sort.multiPhylo()`.

- Update test suite for compatibility with "testthat" > 3.0.4 (@hadley, #83).


# TreeTools 1.5.0

## New functions

- `ConstrainedNJ()` returns an approximation to a neighbour-joining tree
  that respects constraints.

- `PolarizeSplits()` marks a specified taxon as representing the ingroup of all
  splits.

- Add `KeepTip()` and improve performance of `DropTip()`.

- `ImposeConstraint()` makes a tree consistent with topological constraints.

- `as.phylo.Splits()` represents a `Splits` object as a tree.

- `Consensus()` is a faster C++ implementation of `ape::consensus()`.

- `ClusterTable()` C++ functionality imported from "TreeDist".


## Improved functions

- Warn when empty cells passed to `MatrixToPhyDat()`.

- Warn when `LabelSplits(labels)` lack names.

- `SplitFrequency()` drops tips from `forest` that aren't in `reference`.

- `AddTipEverywhere()` supports trees with < 3 leaves.

- Make `RootTree()` and `PhyDatToMatrix()` more robust.

- Support `encoding` option in `ReadCharacters()` function family.

- Support `CHARSTATELABELS` in `ReadCharacters()`.

- Support for more formatting quirks in `ReadNotes()`.

- Better support ambiguous tokens in `WriteTntCharacters()`.


## Optimization

- Fast matching functions from "fastmatch".

- Improve efficiency of `Preorder()` and `Postorder()`,
  and lift limit on tree size.


# TreeTools 1.4.5

- Correct calculation of minimum value in `TCIContext()`.
- Extract tip labels from objects in `StringToPhyDat()`.
- Support `AddTip(tree, where = "tip name")`.
- `SplitFrequency()` supports four-leaf trees.
- Add `RootTree.matrix()` method for edge matrices.
- Add `TipLabels.phyDat()` method.
- Add `NULL` methods for tree reordering functions.
- Additions and improvements to text parsing functions.


# TreeTools 1.4.4

- Add `NTip.phyDat()` method.
- Update `MakeTreeBinary()` docs and tests to reflect updated behaviour of 
  `ape::multi2di()` in 'ape' v5.5.


# TreeTools 1.4.3

 - `AddTip()` supports edge lengths.
 - `CladisticInfo()` supports `Splits` objects.
 - `as.multiPhylo()` converts trees, datasets and Splits objects into 
     `multiPhylo` objects.
 - `LabelSplits(labels = NULL)` labels each split with its associated node.
 - `PhyDatToMatrix()` supports integer-only levels.
 - `SortTree()` supports lists of trees.
 - Improvements to `ReadTntCharacters()` character block extraction
   ([#50](https://github.com/ms609/TreeTools/issues/50)).


# TreeTools 1.4.2

 - Support star trees in `RootTree()`.
 - Improve memory handling in `root_on_node()`.
 - Documentation linkage.


# TreeTools 1.4.1

 - `MSTEdges()` supports distance matrices with > 256 entries.
 - Package 'vdiffr' used conditionally.


# TreeTools 1.4.0

## New functions
 - `MSTLength()` reports length of minimum spanning tree.
 - `AllTipLabels()` returns all labels from all trees in a list.
 - `PairwiseDistances()` (from 'TreeDistData') computes distances between all
   pairs of trees in a list.
 - `ArtificialExtinction()` gains `replaceAll` option.
 - `WriteTntCharacters(types = ...)` writes different character types to TNT 
   file.
 - Tree characterization S3 methods: add `.default` and `.NULL`.

## Enhancements
 - `MSTEdges()` implemented in C++, improving runtime by orders of magnitude.
 - Improved parsing of TNT character files.


# TreeTools 1.3.1

 - Improved parsing of TNT files.
 - Fix misspecified C++ linkage.


# TreeTools 1.3.0

## New functions
 
 - `SisterSize()` and `RootNodeDist()` measure sister-clade size and root-node
   distance.
 - `MSTEdges()`: Edges of minimum spanning tree.
 - `SplitImbalance()`: how balanced is each split?
 - New C++ functions `root_on_node()` and `root_binary()` to root trees quickly
   and robustly.

## Enhancements

 - `TNTReadTree()` handles additional punctuation characters.
 - Import RdMacros package 'Rdpack'.
 
 - C++ implementation of `TipsInSplits()`.
 - Export C++ functions `preorder_edges_and_nodes()` and `postorder_edges()`.
 - Remove obsolete copy of C++ code from 'phangorn'.


# TreeTools 1.2.0

 - `ArtificialExtinction()`: Remove characters that are absent in a fossil 
   template.
 - `WriteTntCharacters()`: Write morphological dataset in TNT format.
 - Improve TNT dataset parsing.
 - Documentation improvements.


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

- Pre-release version spun out of ['TreeSearch'](https://ms609.github.io/TreeSearch/)
  package.  Some functionality is subject to change.
