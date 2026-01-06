# Changelog

## TreeTools 2.0.0.9001 (development)

- Remove hard limit on tree size in `SplitList`.

## TreeTools 2.0.0.9000 (development)

- [`MatrixToPhyDat()`](https://ms609.github.io/TreeTools/dev/reference/MatrixToPhyDat.md)
  gains `tipLabels` parameter.
- Document return value for
  [`J1Index()`](https://ms609.github.io/TreeTools/dev/reference/J1Index.md).

## TreeTools 2.0.0 (2025-09-23)

CRAN release: 2025-09-23

### New functionality

- [`Cherries()`](https://ms609.github.io/TreeTools/dev/reference/Cherries.md)
  counts the cherries in a binary tree.
- New method
  [`as.Splits.integer()`](https://ms609.github.io/TreeTools/dev/reference/Splits.md).
- Add methods for `NULL` objects.

### Fixes and enhancements

- Fix `RoguePlot(sort = TRUE)`
  ([Rogue#33](https://github.com/ms609/Rogue/issues/33)).
- Remove R.cache dependency:
  [`UnrootedKeys()`](https://ms609.github.io/TreeTools/dev/reference/TreeShape.md)
  now uses a native cache implementation.

### Backward incompatible changes

- Require R 3.6.
- Remove deprecated functions `AllDescendantEdges()`,
  `.EnforceOutgroup()`, `ForestSplits()`, `in.Splits()`,
  `PhylogeneticInfo()`, `SpectrumLegend()`, `SplitNumber()`,
  `TreeSplits()`.
- Remove deprecated C++ macro `TREETOOLS_SPLITLIST_INIT`.

## TreeTools 1.16.1 (2025-08-24)

CRAN release: 2025-08-24

- Compiler-safe vector initialization, resolving M1-SAN warnings.

## TreeTools 1.16.0 (2025-08-22)

CRAN release: 2025-08-22

### New functionality

- [`SplitConsistent()`](https://ms609.github.io/TreeTools/dev/reference/SplitConsistent.md)
  calculates split (dis)agreement.
- [`LongBranch()`](https://ms609.github.io/TreeTools/dev/reference/LongBranch.md)
  identifies long-branched taxa.
- [`Treeness()`](https://ms609.github.io/TreeTools/dev/reference/Treeness.md)
  computes the treeness (=stemminess) of a tree, a proxy for its
  phylogenetic signal.
- Add
  [`KeepTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md)
  methods to correspond to
  [`DropTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md).
- [`Preorder()`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md)
  gains `topologyOnly` argument.

### Enhancements

- [`MakeTreeBinary()`](https://ms609.github.io/TreeTools/dev/reference/MakeTreeBinary.md)
  explicitly removes edge lengths.
- Optimize `ClusterTable` class.
- Optimize
  [`NSplits()`](https://ms609.github.io/TreeTools/dev/reference/NSplits.md).
- Remove exported C++ cache objects to prefer calculation with
  intrinsics.
- Other performance improvements.

### Fixes

- Fix
  [`as.ClusterTable()`](https://ms609.github.io/TreeTools/dev/reference/ClusterTable.md)
  when leaf order varies.
- Fix regressions in
  [`as.ClusterTable()`](https://ms609.github.io/TreeTools/dev/reference/ClusterTable.md)
  that caused downstream errors.
- Fix regressions in
  [`PhyToString()`](https://ms609.github.io/TreeTools/dev/reference/PhyToString.md).
- Fix handling of multiple ambiguities in
  [`Reweight()`](https://ms609.github.io/TreeTools/dev/reference/Reweight.md).

## TreeTools 1.15.0 (2025-07-16)

CRAN release: 2025-07-21

- [`Reweight()`](https://ms609.github.io/TreeTools/dev/reference/Reweight.md)
  sets the weight of characters in a phylogenetic dataset.
- [`MatchStrings()`](https://ms609.github.io/TreeTools/dev/reference/MatchStrings.md)
  checks for mismatched tip labels, suggesting corrections.
- [`RootTree()`](https://ms609.github.io/TreeTools/dev/reference/RootTree.md)
  gains `fallback` argument to handle outgroups that do not root a tree.
- Fix
  [`MakeTreeBinary()`](https://ms609.github.io/TreeTools/dev/reference/MakeTreeBinary.md)
  labelling trees as in preorder.
- Fix `as.Splits.matrix(tipLabels != NULL)`.
- Improve performance of
  [`PhyToString()`](https://ms609.github.io/TreeTools/dev/reference/PhyToString.md).
- Modernize aspects of C++ code.

## TreeTools 1.14.0 (2025-05-13)

CRAN release: 2025-05-13

- `AddTip(lengthBelow = NA)` adds leaf at node without adding a new
  edge.

- [`BalancedTree()`](https://ms609.github.io/TreeTools/dev/reference/GenerateTree.md)
  and equivalent gain a `lengths` parameter to specify edge lengths.

- Fix taxa misplaced by `RoguePlot(sort = TRUE)`.

- Fix unexpected polytomies in
  [`Consensus()`](https://ms609.github.io/TreeTools/dev/reference/Consensus.md)
  ([\#168](https://github.com/ms609/TreeTools/issues/168)).

## TreeTools 1.13.1 (2025-04-07)

CRAN release: 2025-04-07

- Support non-unique labels in
  [`DropTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md)
  checks.

## TreeTools 1.13.0 (2025-01-09)

CRAN release: 2025-01-10

- `%in%` and `match` methods for phylo / multiPhylo objects.

- [`Decompose()`](https://ms609.github.io/TreeTools/dev/reference/Decompose.md)
  decomposes additive (ordered) phylogenetic characters by binary
  decomposition.

- Check for overflow in splits functions.

## TreeTools 1.12.0 (2024-07-25)

CRAN release: 2024-07-25

### New methods and functions

- [`TopologyOnly()`](https://ms609.github.io/TreeTools/dev/reference/TopologyOnly.md)
  removes metadata from phylo objects.
- [`J1Index()`](https://ms609.github.io/TreeTools/dev/reference/J1Index.md)
  computes the robust, universal tree balance measure of Lemant *et al*.
  2022 <doi:10.1093/sysbio/syac027>, incorporating
  [code](https://github.com/robjohnnoble/RUtreebalance) by Rob Noble.

### Enhancements

- Consistent sequence of list entries in phylo objects.
- [`RandomTree()`](https://ms609.github.io/TreeTools/dev/reference/GenerateTree.md)
  returns trees for \< 3 leaves.
- [`root_on_node()`](https://ms609.github.io/TreeTools/dev/reference/root_on_node.md)
  handles trees with \< 2 leaves.
- Support larger trees in
  [`TotalCopheneticIndex()`](https://ms609.github.io/TreeTools/dev/reference/TotalCopheneticIndex.md),
  fixing [\#158](https://github.com/ms609/TreeTools/issues/158).

## TreeTools 1.11.1 (2024-06-06)

CRAN release: 2024-06-07

- Set random seed and increase tolerance to avoid false negatives on
  tests.

## TreeTools 1.11.0 (2024-05-23)

CRAN release: 2024-06-05

### New methods and functions

- [`YuleTree()`](https://ms609.github.io/TreeTools/dev/reference/GenerateTree.md)
  generates a random tree by the Yule process.
- [`DescendantTips()`](https://ms609.github.io/TreeTools/dev/reference/DescendantEdges.md)
  complements
  [`DescendantEdges()`](https://ms609.github.io/TreeTools/dev/reference/DescendantEdges.md),
  rewritten in C++, fixing a bug when edges were not in preorder.
- [`NodeNumbers()`](https://ms609.github.io/TreeTools/dev/reference/NodeNumbers.md)
  returns the indices of nodes within a tree.

### Enhancements

- `RandomTree(root = TRUE)` roots the tree on a random edge.
- `RoguePlot()$legendLabels` returns suggested labels for legend.
- Support node labels in
  [`AddTip()`](https://ms609.github.io/TreeTools/dev/reference/AddTip.md),
  [`CollapseNode()`](https://ms609.github.io/TreeTools/dev/reference/CollapseNode.md),
  [`DropTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md),
  [`MakeTreeBinary()`](https://ms609.github.io/TreeTools/dev/reference/MakeTreeBinary.md),
  [`Renumber()`](https://ms609.github.io/TreeTools/dev/reference/Renumber.md),
  [`Reorder()`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md),
  [`SortTree()`](https://ms609.github.io/TreeTools/dev/reference/SortTree.md),
  [`Subtree()`](https://ms609.github.io/TreeTools/dev/reference/Subtree.md)
  ([\#149](https://github.com/ms609/TreeTools/issues/149)).
- `AddTip(edgeLength = NULL)` defaults to `lengthBelow`. This will
  become the default in a future release.
- An entry point to the C++ function
  [`root_on_node()`](https://ms609.github.io/TreeTools/dev/reference/root_on_node.md)
  is now exported (intended for expert use only).
- Fix handling of weighted trees by
  [`root_on_node()`](https://ms609.github.io/TreeTools/dev/reference/root_on_node.md).
- Use
  [`KeepTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md)
  internally so
  [`SplitFrequency()`](https://ms609.github.io/TreeTools/dev/reference/SplitFrequency.md)
  supports `Splits` objects as documented.

## TreeTools 1.10.0 (2023-08-18)

CRAN release: 2023-08-18

### New methods and functions

- [`TipTimedTree()`](https://ms609.github.io/TreeTools/dev/reference/TipTimedTree.md)
  displays trees where leaves are associated with absolute ages.
- [`ReadMrBayesTrees()`](https://ms609.github.io/TreeTools/dev/reference/ReadMrBayesTrees.md)
  samples trees from posterior of MrBayes output.
- [`is.TreeNumber()`](https://ms609.github.io/TreeTools/dev/reference/is.TreeNumber.md)
  method.

### Improvements

- Support zero-edge trees in
  [`as.Splits()`](https://ms609.github.io/TreeTools/dev/reference/Splits.md)
  and
  [`NSplits()`](https://ms609.github.io/TreeTools/dev/reference/NSplits.md).
- Support empty constraints in
  [`AddUnconstrained()`](https://ms609.github.io/TreeTools/dev/reference/ImposeConstraint.md).
- Add space between tokens in
  [`WriteTntCharacters()`](https://ms609.github.io/TreeTools/dev/reference/WriteTntCharacters.md)
  to support continuous characters
  ([\#139](https://github.com/ms609/TreeTools/issues/139)).

### Deprecations and breaking changes

- Change order of parameters in
  [`DescendantEdges()`](https://ms609.github.io/TreeTools/dev/reference/DescendantEdges.md)
- Deprecate `AllDescendantEdges()`; use
  [`DescendantEdges()`](https://ms609.github.io/TreeTools/dev/reference/DescendantEdges.md)
  instead.
- Deprecate `EnforceOutgroup()`; use
  [`RootTree()`](https://ms609.github.io/TreeTools/dev/reference/RootTree.md)
  instead.
- Remove `NonDuplicateRoot()` and `in.Splits()`.

## TreeTools 1.9.2 (2023-04-25)

CRAN release: 2023-04-27

- Improve support for comments in
  [`ReadNotes()`](https://ms609.github.io/TreeTools/dev/reference/ReadCharacters.md).
- Support Nexus-escaped `''`s in
  [`ReadCharacters()`](https://ms609.github.io/TreeTools/dev/reference/ReadCharacters.md).
- Add `legend` parameter to
  [`RoguePlot()`](https://ms609.github.io/TreeTools/dev/reference/RoguePlot.md).
- [`RoguePlot()`](https://ms609.github.io/TreeTools/dev/reference/RoguePlot.md)
  now returns invisibly.
- Deprecate `SpectrumLegend()` – spun off to separate
  [“PlotTools”](https://ms609.github.io/PlotTools/) package.

## TreeTools 1.9.1 (2023-03-21)

CRAN release: 2023-03-20

- [`AddUnconstrained()`](https://ms609.github.io/TreeTools/dev/reference/ImposeConstraint.md)
  and
  [`ImposeConstraint()`](https://ms609.github.io/TreeTools/dev/reference/ImposeConstraint.md)
  handle wider range of inputs.

- [`PhyDatToMatrix()`](https://ms609.github.io/TreeTools/dev/reference/MatrixToPhyDat.md)
  can (and by default does) override levels to write ambiguous tokens in
  custom formats such as `{01}`.

- Call C functions using symbols, not strings.

## TreeTools 1.9.0 (2022-11-29)

CRAN release: 2022-11-28

### New methods and functions

- [`ZeroTaxonTree()`](https://ms609.github.io/TreeTools/dev/reference/TrivialTree.md)
  creates a `phylo` object with no leaves.

- [`DropTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md)
  gains new methods
  [`DropTip.list()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md)
  and
  [`DropTip.NULL()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md).

- `as.matrix.phylo()` converts a tree to a matrix representation,
  allowing a tree to be passed as a constraint to
  [`ImposeConstraint()`](https://ms609.github.io/TreeTools/dev/reference/ImposeConstraint.md).

- `as.matrix.Splits()` and `as.matrix.phyDat()` methods added as
  synonyms to
  [`as.logical.Splits()`](https://ms609.github.io/TreeTools/dev/reference/Splits.md)
  and
  [`PhyDatToMatrix()`](https://ms609.github.io/TreeTools/dev/reference/MatrixToPhyDat.md).

### Improvements

- Handle `TipLabels(0)` and `BalancedTree(0)`.

- Support zero-leaf trees in
  [`as.Splits()`](https://ms609.github.io/TreeTools/dev/reference/Splits.md)
  and `duplicated.Splits()`.

- Support non-identical tip labels in
  [`as.Splits()`](https://ms609.github.io/TreeTools/dev/reference/Splits.md).

- Try Latin-1 encoding if
  [`ReadCharacters()`](https://ms609.github.io/TreeTools/dev/reference/ReadCharacters.md)
  family fail under UTF-8.

## TreeTools 1.8.0 (2022-09-15)

CRAN release: 2022-09-15

### New methods and functions

- [`TntOrder()`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md)
  renumbers a tree’s nodes to match TNT’s convention.

- [`head()`](https://rdrr.io/r/utils/head.html) and
  [`tail()`](https://rdrr.io/r/utils/head.html) methods for Splits
  objects.

- Set names of splits object with `names(splits) <- ...`.

- [`as.Splits()`](https://ms609.github.io/TreeTools/dev/reference/Splits.md)
  support character vectors in the form “…\*\*\*“.

### Improvements

- [`ReadTntTree()`](https://ms609.github.io/TreeTools/dev/reference/ReadTntTree.md)
  reads tree tags and follows TNT node numbering conventions.

- `SpectrumLegend()` gains `title` parameter and more styling options.

- Support \> 32767 trees in
  [`Consensus()`](https://ms609.github.io/TreeTools/dev/reference/Consensus.md)
  ([\#127](https://github.com/ms609/TreeTools/issues/127)).

- [`DropTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md)
  speed improved when branch lengths are present.

## TreeTools 1.7.3 (2022-07-20)

CRAN release: 2022-07-20

- [`ReadTntTree()`](https://ms609.github.io/TreeTools/dev/reference/ReadTntTree.md)
  supports multi-line trees.

- [`as.MixedBase()`](https://ms609.github.io/TreeTools/dev/reference/TreeNumber.md)
  supports larger trees (44-32767 tips).

- Add deprecation warning to `in.Splits()`.

## TreeTools 1.7.2 (2022-05-24)

CRAN release: 2022-05-24

- [`RenumberTips()`](https://ms609.github.io/TreeTools/dev/reference/RenumberTips.md)
  drops “preorder” attribute, as reordering tip labels may break edge
  ordering guarantee.

- Native implementation of `ClusterTable` class.

- Replace `throw` with `stop` in C++ scripts.

## TreeTools 1.7.1 (2022-03-25)

CRAN release: 2022-03-25

- [`AddTip()`](https://ms609.github.io/TreeTools/dev/reference/AddTip.md):
  Fix bug when adding tip to root of weighted tree.

## TreeTools 1.7.0 (2022-03-21)

CRAN release: 2022-03-21

### New methods and functions

- `rev.Splits()` reverses order in which splits are listed.

- [`KeepTip.Splits()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md)
  is a faster alternative to `SubSplit()`.

- `%in%.Splits()` retains names when comparing small splits
  ([\#40](https://github.com/ms609/TreeTools/issues/40)).

- [`sort.multiPhylo()`](https://ms609.github.io/TreeTools/dev/reference/sort.multiPhylo.md)
  sorts lists of trees according to their mixed base representation
  ([\#84](https://github.com/ms609/TreeTools/issues/84)).

- Bitwise manipulation of splits with `|`, `&`, `xor`.

- [`as.MixedBase()`](https://ms609.github.io/TreeTools/dev/reference/TreeNumber.md)
  uniquely represents binary trees as a mixed-base vector.

- [`PathLengths()`](https://ms609.github.io/TreeTools/dev/reference/PathLengths.md)
  describes all paths within a tree.

- [`KeptVerts()`](https://ms609.github.io/TreeTools/dev/reference/KeptVerts.md)
  and
  [`KeptPaths()`](https://ms609.github.io/TreeTools/dev/reference/KeptPaths.md)
  identify elements in reduced trees.

- [`PostorderOrder()`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md)
  describes a sequence of edges corresponding to a postorder traversal
  of a tree.

- `SpectrumLegend()` adds gradients to plot legends.

### Improvements

- Improve handling of zero-split trees.

- [`DropTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md)
  no longer adds a root to unrooted trees, and retains edge lengths.

- Improve speed of
  [`DropTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md),
  by an order of magnitude in some cases.

- Support edge lengths in
  [`Preorder()`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md),
  [`RootTree()`](https://ms609.github.io/TreeTools/dev/reference/RootTree.md),
  [`UnrootTree()`](https://ms609.github.io/TreeTools/dev/reference/RootTree.md)
  and
  [`Postorder()`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md)
  ([\#49](https://github.com/ms609/TreeTools/issues/49),
  [\#89](https://github.com/ms609/TreeTools/issues/89)).

- Fix bug when tree is rooted on a discontinuous outgroup.

- [`SortTree()`](https://ms609.github.io/TreeTools/dev/reference/SortTree.md)
  handles weighted and non-binary trees
  ([\#25](https://github.com/ms609/TreeTools/issues/25),
  [\#25](https://github.com/ms609/TreeTools/issues/49)), and gains
  option to sort by tip labels.

- `TipsInSplits(smallest = TRUE)` counts tips in smaller bipartition.

- Fix a bug with `phyDat` objects in
  [`ArtificialExtinction()`](https://ms609.github.io/TreeTools/dev/reference/ArtificialExtinction.md).

- [`RenumberTips()`](https://ms609.github.io/TreeTools/dev/reference/RenumberTips.md)
  allows `tipOrder` to contain elements not present in `tree`.

- Use lighter Rcpp headers.

- Small improvements to computational efficiency.

### Deprecations

- Remove deprecated function `PostorderEdges()`
  ([\#35](https://github.com/ms609/TreeTools/issues/35)).

## TreeTools 1.6.0 (2022-01-12)

CRAN release: 2022-01-12

### New functions

- [`RoguePlot()`](https://ms609.github.io/TreeTools/dev/reference/RoguePlot.md)
  plots the positions of rogue taxa.

### Improvements

- [`DropTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md)
  gains `check` parameter to allow slightly faster operation where input
  is guaranteed to be valid.

- [`RandomTree()`](https://ms609.github.io/TreeTools/dev/reference/GenerateTree.md)
  gains `nodes` parameter allow the inclusion of polytomies.

- Infer `tips` parameter if missing in
  [`StringToPhyDat()`](https://ms609.github.io/TreeTools/dev/reference/PhyToString.md).

- Remove dependency on “phangorn” (allowing use on R \< 4.1)

- Improve parsing of information from nexus files.

- Export
  [`DropTipPhylo()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md)
  as wrapper to
  [`DropTip.phylo()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md).

## TreeTools 1.5.1 (2021-10-06)

CRAN release: 2021-10-06

- [`PhyDatToMatrix()`](https://ms609.github.io/TreeTools/dev/reference/MatrixToPhyDat.md)
  optionally encodes ambiguous / inapplicable tokens as `NA`.

- Implement
  [`sort.multiPhylo()`](https://ms609.github.io/TreeTools/dev/reference/sort.multiPhylo.md).

- Update test suite for compatibility with “testthat” \> 3.0.4
  ([@hadley](https://github.com/hadley),
  [\#83](https://github.com/ms609/TreeTools/issues/83)).

## TreeTools 1.5.0 (2021-09-13)

CRAN release: 2021-09-08

### New functions

- [`ConstrainedNJ()`](https://ms609.github.io/TreeTools/dev/reference/ConstrainedNJ.md)
  returns an approximation to a neighbour-joining tree that respects
  constraints.

- [`PolarizeSplits()`](https://ms609.github.io/TreeTools/dev/reference/PolarizeSplits.md)
  marks a specified taxon as representing the ingroup of all splits.

- Add
  [`KeepTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md)
  and improve performance of
  [`DropTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md).

- [`ImposeConstraint()`](https://ms609.github.io/TreeTools/dev/reference/ImposeConstraint.md)
  makes a tree consistent with topological constraints.

- `as.phylo.Splits()` represents a `Splits` object as a tree.

- [`Consensus()`](https://ms609.github.io/TreeTools/dev/reference/Consensus.md)
  is a faster C++ implementation of
  [`ape::consensus()`](https://rdrr.io/pkg/ape/man/consensus.html).

- [`ClusterTable()`](https://ms609.github.io/TreeTools/dev/reference/ClusterTable.md)
  C++ functionality imported from “TreeDist”.

### Improved functions

- Warn when empty cells passed to
  [`MatrixToPhyDat()`](https://ms609.github.io/TreeTools/dev/reference/MatrixToPhyDat.md).

- Warn when `LabelSplits(labels)` lack names.

- [`SplitFrequency()`](https://ms609.github.io/TreeTools/dev/reference/SplitFrequency.md)
  drops tips from `forest` that aren’t in `reference`.

- [`AddTipEverywhere()`](https://ms609.github.io/TreeTools/dev/reference/AddTip.md)
  supports trees with \< 3 leaves.

- Make
  [`RootTree()`](https://ms609.github.io/TreeTools/dev/reference/RootTree.md)
  and
  [`PhyDatToMatrix()`](https://ms609.github.io/TreeTools/dev/reference/MatrixToPhyDat.md)
  more robust.

- Support `encoding` option in
  [`ReadCharacters()`](https://ms609.github.io/TreeTools/dev/reference/ReadCharacters.md)
  function family.

- Support `CHARSTATELABELS` in
  [`ReadCharacters()`](https://ms609.github.io/TreeTools/dev/reference/ReadCharacters.md).

- Support for more formatting quirks in
  [`ReadNotes()`](https://ms609.github.io/TreeTools/dev/reference/ReadCharacters.md).

- Better support ambiguous tokens in
  [`WriteTntCharacters()`](https://ms609.github.io/TreeTools/dev/reference/WriteTntCharacters.md).

### Optimization

- Fast matching functions from “fastmatch”.

- Improve efficiency of
  [`Preorder()`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md)
  and
  [`Postorder()`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md),
  and lift limit on tree size.

## TreeTools 1.4.5 (2021-06-23)

CRAN release: 2021-06-23

- Correct calculation of minimum value in
  [`TCIContext()`](https://ms609.github.io/TreeTools/dev/reference/TotalCopheneticIndex.md).
- Extract tip labels from objects in
  [`StringToPhyDat()`](https://ms609.github.io/TreeTools/dev/reference/PhyToString.md).
- Support `AddTip(tree, where = "tip name")`.
- [`SplitFrequency()`](https://ms609.github.io/TreeTools/dev/reference/SplitFrequency.md)
  supports four-leaf trees.
- Add `RootTree.matrix()` method for edge matrices.
- Add
  [`TipLabels.phyDat()`](https://ms609.github.io/TreeTools/dev/reference/TipLabels.md)
  method.
- Add `NULL` methods for tree reordering functions.
- Additions and improvements to text parsing functions.

## TreeTools 1.4.4 (2021-04-23)

CRAN release: 2021-04-20

- Add
  [`NTip.phyDat()`](https://ms609.github.io/TreeTools/dev/reference/NTip.md)
  method.
- Update
  [`MakeTreeBinary()`](https://ms609.github.io/TreeTools/dev/reference/MakeTreeBinary.md)
  docs and tests to reflect updated behaviour of
  [`ape::multi2di()`](https://rdrr.io/pkg/ape/man/multi2di.html) in
  ‘ape’ v5.5.

## TreeTools 1.4.3 (2021-04-12)

CRAN release: 2021-04-12

- [`AddTip()`](https://ms609.github.io/TreeTools/dev/reference/AddTip.md)
  supports edge lengths.
- [`CladisticInfo()`](https://ms609.github.io/TreeTools/dev/reference/CladisticInfo.md)
  supports `Splits` objects.
- [`as.multiPhylo()`](https://ms609.github.io/TreeTools/dev/reference/as.multiPhylo.md)
  converts trees, datasets and Splits objects into `multiPhylo` objects.
- `LabelSplits(labels = NULL)` labels each split with its associated
  node.
- [`PhyDatToMatrix()`](https://ms609.github.io/TreeTools/dev/reference/MatrixToPhyDat.md)
  supports integer-only levels.
- [`SortTree()`](https://ms609.github.io/TreeTools/dev/reference/SortTree.md)
  supports lists of trees.
- Improvements to
  [`ReadTntCharacters()`](https://ms609.github.io/TreeTools/dev/reference/ReadCharacters.md)
  character block extraction
  ([\#50](https://github.com/ms609/TreeTools/issues/50)).

## TreeTools 1.4.2 (2021-01-26)

CRAN release: 2021-01-26

- Support star trees in
  [`RootTree()`](https://ms609.github.io/TreeTools/dev/reference/RootTree.md).
- Improve memory handling in
  [`root_on_node()`](https://ms609.github.io/TreeTools/dev/reference/root_on_node.md).
- Documentation linkage.

## TreeTools 1.4.1 (2020-12-09)

CRAN release: 2020-12-09

- [`MSTEdges()`](https://ms609.github.io/TreeTools/dev/reference/MSTEdges.md)
  supports distance matrices with \> 256 entries.
- Package ‘vdiffr’ used conditionally.

## TreeTools 1.4.0 (2020-10-20)

CRAN release: 2020-10-19

### New functions

- [`MSTLength()`](https://ms609.github.io/TreeTools/dev/reference/MSTEdges.md)
  reports length of minimum spanning tree.
- [`AllTipLabels()`](https://ms609.github.io/TreeTools/dev/reference/TipLabels.md)
  returns all labels from all trees in a list.
- [`PairwiseDistances()`](https://ms609.github.io/TreeTools/dev/reference/PairwiseDistances.md)
  (from ‘TreeDistData’) computes distances between all pairs of trees in
  a list.
- [`ArtificialExtinction()`](https://ms609.github.io/TreeTools/dev/reference/ArtificialExtinction.md)
  gains `replaceAll` option.
- `WriteTntCharacters(types = ...)` writes different character types to
  TNT file.
- Tree characterization S3 methods: add `.default` and `.NULL`.

### Enhancements

- [`MSTEdges()`](https://ms609.github.io/TreeTools/dev/reference/MSTEdges.md)
  implemented in C++, improving runtime by orders of magnitude.
- Improved parsing of TNT character files.

## TreeTools 1.3.1

CRAN release: 2020-10-03

- Improved parsing of TNT files.
- Fix misspecified C++ linkage.

## TreeTools 1.3.0 (2020-09-22)

CRAN release: 2020-09-22

### New functions

- [`SisterSize()`](https://ms609.github.io/TreeTools/dev/reference/Stemwardness.md)
  and
  [`RootNodeDist()`](https://ms609.github.io/TreeTools/dev/reference/Stemwardness.md)
  measure sister-clade size and root-node distance.
- [`MSTEdges()`](https://ms609.github.io/TreeTools/dev/reference/MSTEdges.md):
  Edges of minimum spanning tree.
- [`SplitImbalance()`](https://ms609.github.io/TreeTools/dev/reference/TipsInSplits.md):
  how balanced is each split?
- New C++ functions
  [`root_on_node()`](https://ms609.github.io/TreeTools/dev/reference/root_on_node.md)
  and `root_binary()` to root trees quickly and robustly.

### Enhancements

- `TNTReadTree()` handles additional punctuation characters.

- Import RdMacros package ‘Rdpack’.

- C++ implementation of
  [`TipsInSplits()`](https://ms609.github.io/TreeTools/dev/reference/TipsInSplits.md).

- Export C++ functions `preorder_edges_and_nodes()` and
  `postorder_edges()`.

- Remove obsolete copy of C++ code from ‘phangorn’.

## TreeTools 1.2.0 (2020-08-30)

CRAN release: 2020-08-03

- [`ArtificialExtinction()`](https://ms609.github.io/TreeTools/dev/reference/ArtificialExtinction.md):
  Remove characters that are absent in a fossil template.
- [`WriteTntCharacters()`](https://ms609.github.io/TreeTools/dev/reference/WriteTntCharacters.md):
  Write morphological dataset in TNT format.
- Improve TNT dataset parsing.
- Documentation improvements.

## TreeTools 1.1.0 (2020-07-07)

CRAN release: 2020-07-07

### New functions

- [`RandomTree()`](https://ms609.github.io/TreeTools/dev/reference/GenerateTree.md):
  Draw tree from uniform distribution, instead of via
  [`ape::rtree()`](https://rdrr.io/pkg/ape/man/rtree.html).
- [`MakeTreeBinary()`](https://ms609.github.io/TreeTools/dev/reference/MakeTreeBinary.md):
  Uniform equivalent of
  [`ape::multi2di()`](https://rdrr.io/pkg/ape/man/multi2di.html).
- `match.list()` method for lists of splits.
- [`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/dev/reference/SplitsInBinaryTree.md):
  How many splits occur in an *n*-leaf binary tree?
- [`vapply64()`](https://ms609.github.io/TreeTools/dev/reference/sapply64.md),
  [`sapply64()`](https://ms609.github.io/TreeTools/dev/reference/sapply64.md),
  [`replicate64()`](https://ms609.github.io/TreeTools/dev/reference/sapply64.md):
  helper functions when a function returns a 64-bit integer.

### Enhancements

- Use methods for
  [`UnrootTree()`](https://ms609.github.io/TreeTools/dev/reference/RootTree.md),
  [`RootTree()`](https://ms609.github.io/TreeTools/dev/reference/RootTree.md),
  [`RootOnNode()`](https://ms609.github.io/TreeTools/dev/reference/RootTree.md)
  to support lists of trees.

## TreeTools 1.0.0 (2020-06-08)

CRAN release: 2020-06-08

### New functions

- [`CladisticInfo()`](https://ms609.github.io/TreeTools/dev/reference/CladisticInfo.md):
  Calculate the information content of a tree.
- [`RootNode()`](https://ms609.github.io/TreeTools/dev/reference/RootNode.md):
  Which node is a tree’s root?
- [`UnrootTree()`](https://ms609.github.io/TreeTools/dev/reference/RootTree.md):
  Safely remove a root node.
- [`NodeDepth()`](https://ms609.github.io/TreeTools/dev/reference/NodeDepth.md):
  Discriminate shallow from deep nodes.
- [`NodeOrder()`](https://ms609.github.io/TreeTools/dev/reference/NodeOrder.md),
  [`NDescendants()`](https://ms609.github.io/TreeTools/dev/reference/NDescendants.md):
  Count edges incident to each node.
- [`CladeSizes()`](https://ms609.github.io/TreeTools/dev/reference/CladeSizes.md):
  Count leaves / nodes descended from each node.
- [`ListAncestors()`](https://ms609.github.io/TreeTools/dev/reference/ListAncestors.md):
  List ancestors of a node.
- [`LabelSplits()`](https://ms609.github.io/TreeTools/dev/reference/LabelSplits.md):
  Label splits on plotted tree.
- [`DropTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md):
  Remove tip, handling weird node orders.
- [`LeafLabelInterchange()`](https://ms609.github.io/TreeTools/dev/reference/LeafLabelInterchange.md):
  Exchange position of *n* tips.
- [`StarTree()`](https://ms609.github.io/TreeTools/dev/reference/GenerateTree.md):
  Generate unresolved tree.
- [`TotalCopheneticIndex()`](https://ms609.github.io/TreeTools/dev/reference/TotalCopheneticIndex.md)
  integrated from ‘tci’ package.

### Deprecations

- `PostorderEdges()`: use
  [`Postorder()`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md)
  instead.
- `NonDuplicateRoot()`: unused internal function.
- `match.Splits()`: use
  [`match()`](https://ms609.github.io/TreeTools/dev/reference/match.Splits.md)
  instead.
- `in.Splits()`: use `%in%.Splits()` instead.

### Enhancements

- Improve support for unrooted trees in
  [`as.Splits()`](https://ms609.github.io/TreeTools/dev/reference/Splits.md).
- Use methods so `Reorder` functions can handle `multiPhylo` objects and
  edges.
- Handle funny node orders.
- Support continuous characters in
  [`ReadCharacters()`](https://ms609.github.io/TreeTools/dev/reference/ReadCharacters.md).
- Improve performance of
  [`as.logical.Splits()`](https://ms609.github.io/TreeTools/dev/reference/Splits.md)
  and related functions.
- Fail nicely when trees are too large for memory.
- Fix memory leak in
  [`as.Splits()`](https://ms609.github.io/TreeTools/dev/reference/Splits.md).
- Various under-the-hood improvements to functions.
- Documentation improvements.

## TreeTools 0.1.4 (2020-03-04)

CRAN release: 2020-03-04

- Catch hang-inducing bugs in
  [`RootOnNode()`](https://ms609.github.io/TreeTools/dev/reference/RootTree.md).
- Update `doubleFactorials` cache to fix
  [`as.integer()`](https://rdrr.io/r/base/integer.html) rounding error.
- Support unrooted trees in
  [`AddTipEverywhere()`](https://ms609.github.io/TreeTools/dev/reference/AddTip.md).
- Documentation improvements.

## TreeTools 0.1.3 (2020-01-07)

CRAN release: 2019-12-19

- [`RootOnNode()`](https://ms609.github.io/TreeTools/dev/reference/RootTree.md):
  Quickly root a tree on a specified node.
- Improve portability of C++ code.

## TreeTools 0.1.2 (2020-12-18)

CRAN release: 2019-12-18

- `as.Newick`: Fast conversion to Newick format.
- `as.TreeNumber`: Tree shape enumeration.

## TreeTools 0.1.1

- Add functions to translate trees to mixed base integers.
- `RenumberTips` can extract tip order from `phylo` and `Splits`
  objects.
- Documentation changes to satisfy CRAN submission requirements.

## TreeTools 0.1.0 (2019-10-30)

- Pre-release version spun out of
  [‘TreeSearch’](https://ms609.github.io/TreeSearch/) package. Some
  functionality is subject to change.
