#' Number of trees matching a bipartition split
#'
#' Calculates the number of unrooted bifurcated trees that are consistent with
#' a bipartition split that divides taxa into groups of size `A` and `B`.
#'
#' @param A,B Integer specifying the number of taxa in each partition.
#'
#' @return `TreesMatchingSplit()` returns a numeric specifying the number of trees
#' that are compatible with the given split.
#'
#' `LnTreesMatchingSplit()` and `Log2TreesMatchingSplit()` give the natural
#' and base-2 logarithms of this number.
#'
#' @examples
#' TreesMatchingSplit(5, 6)
#' LnTreesMatchingSplit(5, 6)
#' Log2TreesMatchingSplit(5, 6)
#'
#' @template MRS
#'
#' @family split information functions
#' @export
TreesMatchingSplit <- function(A, B = A[2]) {
  if (A[1] == 0) NUnrooted(B) else
  if (B == 0) NUnrooted(A[1]) else
  NRooted(A[1]) * NRooted(B)
}

#' @rdname TreesMatchingSplit
#' @export
LnTreesMatchingSplit <- function(A, B = A[2]) {
  if (A[1] == 0) LnUnrooted.int(B) else
  if (B == 0) LnUnrooted.int(A[1]) else
  LnRooted.int(A[1]) + LnRooted.int(B)
}

#' @rdname TreesMatchingSplit
#' @export
Log2TreesMatchingSplit <- function(A, B = A[2]) {
  if (A[1] == 0) Log2Unrooted.int(B) else
  if (B == 0) Log2Unrooted.int(A[1]) else
  Log2Rooted.int(A[1]) + Log2Rooted.int(B)
}

#' Character information content
#'
#' `CharacterInformation()` calculates the cladistic information content
#' \insertCite{Steel2006}{TreeTools} of a given character, in bits.
#' The total information in all characters gives a measure of the potential
#' utility of a dataset \insertCite{Cotton2008}{TreeTools}, which can be
#' compared with a profile parsimony score \insertCite{Faith2001}{TreeTools} to
#' evaluate the degree of homoplasy within a dataset.
#'
#' @param tokens Character vector specifying the tokens assigned to each taxon for
#' a character.  Example: `c(0, 0, 0, 1, 1, 1, '?', '-')`.
#'
#' Note that ambiguous tokens such as `(01)` are not supported, and should be
#' replaced with `?`.
#'
#' @return `CharacterInformation()` returns a numeric specifying the
#' phylogenetic information content of the character (_sensu_ 
#' \insertCite{Steel2006;nobrackets}{TreeTools}), in bits.
#'
#' @references
#' - \insertAllCited{}
#' @family split information functions
#' @template MRS
#' @importFrom fastmatch %fin%
#' @export
CharacterInformation <- function(tokens) {
  tokenCounts <- table(tokens)
  # Our character splits our taxa into groups with the same token
  # ?s and -s are best ignored
  splits <- tokenCounts[!(names(tokenCounts) %fin% c('?', '-'))]

  # Return:
  MultiSplitInformation(splits)
}

#' Phylogenetic information content of splitting leaves into two partitions
#'
#' Calculate the phylogenetic information content (_sensu_ 
#' \insertCite{Steel2006;nobrackets}{TreeTools}) of a split, which
#' reflects the probability that a uniformly selected random tree will contain#
#' the split: a split that is consistent with a smaller number of trees will
#' have a higher information content.
#'
#' `SplitInformation()` addresses bipartition splits, which correspond to
#' edges in an unrooted phylogeny; `MultiSplitInformation()` supports splits
#' that subdivide taxa into multiple partitions, which may correspond to
#' multi-state characters in a phylogenetic matrix.
#'
#' A simple way to characterise trees is to count the number of edges.
#' (Edges are almost, but not quite, equivalent to nodes.)
#' Counting edges (or nodes) provides a quick measure of a tree's resolution,
#' and underpins the Robinson-Foulds tree distance measure.
#' Not all edges, however, are created equal.
#'
#' An edge splits the leaves of a tree into two subdivisions.  The more equal
#' these subdivisions are in size, the more instructive this edge is.
#' Intuitively, the division of mammals from reptiles is a profound revelation
#' that underpins much of zoology; recognizing that two species of bat are more
#' closely related to each other than to any other mammal or reptile is still
#' instructive, but somewhat less fundamental.
#'
#' Formally, the phylogenetic (Shannon) information content of a split _S_,
#' _h(S)_, corresponds to the probability that a uniformly selected random tree
#' will contain the split, _P(S)_: _h(S)_ = -log _P(S)_.
#' Base 2 logarithms are typically employed to yield an information content in
#' bits.
#'
#' As an example, the split `AB|CDEF` occurs in 15 of the 105 six-leaf trees;
#' _h_(`AB|CDEF`) = -log _P_(`AB|CDEF`) = -log(15/105) ~ 2.81 bits.  The split
#' `ABC|DEF` subdivides the leaves more evenly, and is thus more instructive:
#' it occurs in just nine of the 105 six-leaf trees, and
#' _h_(`ABC|DEF`) = -log(9/105) ~ 3.54 bits.
#'
#' As the number of leaves increases, a single even split may contain more
#' information than multiple uneven splits -- see the examples section below.
#'
#' Summing the information content of all splits within a tree, perhaps using
#' the '[TreeDist](https://ms609.github.io/TreeDist/)' function
#' [`SplitwiseInfo()`](https://ms609.github.io/TreeDist/reference/TreeInfo.html),
#' arguably gives a more instructive picture of its resolution than simply
#' counting the number of splits that are present -- though with the caveat
#' that splits within a tree are not independent of one another, so some
#' information may be double counted.  (This same charge applies to simply
#' counting nodes, too.)
#'
#' Alternatives would be to count the number of quartets that are resolved,
#' perhaps using the '[Quartet](https://ms609.github.io/Quartet/)' function
#' [`QuartetStates()`](https://ms609.github.io/Quartet/reference/QuartetState.html),
#' or to use a different take on the information contained within a split, the
#' clustering information: see the 'TreeDist' function
#' [`ClusteringInfo()`](https://ms609.github.io/TreeDist/reference/TreeInfo.html)
#' for details.
#'
#'
#' @inheritParams TreesMatchingSplit
#'
#' @return `SplitInformation()` and `MultiSplitInformation()` return the
#' phylogenetic information content, in bits, of a split that subdivides leaves
#' into partitions of the specified sizes.
#'
#' @examples
#' # Eight leaves can be split evenly:
#' SplitInformation(4, 4)
#'
#' # or unevenly, which is less informative:
#' SplitInformation(2, 6)
#'
#' # A single split that evenly subdivides 50 leaves contains more information
#' # that seven maximally uneven splits on the same leaves:
#' SplitInformation(25, 25)
#' 7 * SplitInformation(2, 48)
#' @references \insertAllCited{}
#'
#' @family split information functions
#'
#' @seealso
#'
#' Sum the phylogenetic information content of splits within a tree:
#'  [`TreeDist::SplitwiseInfo()`](https://ms609.github.io/TreeDist/reference/TreeInfo.html)
#'
#'
#' Sum the clustering information content of splits within a tree:
#'  [`TreeDist::ClusteringInfo()`](https://ms609.github.io/TreeDist/reference/TreeInfo.html)
#'
#'
#'
#' @template MRS
#' @export
SplitInformation <- function(A, B = A[1]) {
  -(Log2TreesMatchingSplit(A, B) - Log2Unrooted.int(A + B))
}

#' @rdname SplitInformation
#' @param partitionSizes Integer vector specifying the number of taxa in each
#' partition of a multi-partition split.
#' @examples
#' # Three ways to split eight leaves into multiple partitions:
#' MultiSplitInformation(c(2, 2, 4))
#' MultiSplitInformation(c(2, 3, 3))
#' MultiSplitInformation(rep(2, 4))
#'
#'
#' @export
MultiSplitInformation <- function(partitionSizes) {
  Log2Unrooted.int(sum(partitionSizes)) - Log2UnrootedMult(partitionSizes)
}

#' Number of trees consistent with split
#'
#' Calculates the number of unrooted bifurcating trees consistent with the
#' specified multi-partition split, using theorem two of
#' \insertCite{Carter1990;textual}{TreeTools}.
#'
#' @param \dots A series or vector of integers listing the number of tips in
#' each of a number of tree splits (e.g. bipartitions).
#' For example, `3, 5` states that a character divides a set of eight tips into
#' a group of three and a group of five.
#TODO Support Splits objects, using a Method.
#'
#' @return `UnrootedTreesMatchingSplit()` returns an integer specifying the
#'  number of unrooted bifurcating trees consistent with the specified split.
#'
#'
#' @references \insertAllCited{}
#' 
#' @examples
#' UnrootedTreesMatchingSplit(c(3, 5))
#' UnrootedTreesMatchingSplit(3, 2, 1, 2)
#' @template MRS
#' @family split information functions
#' @export
UnrootedTreesMatchingSplit <- function(...) {
  # use exp and log as it's just as fast, but less likely to overflow to Inf
  exp(LnUnrootedTreesMatchingSplit(...))
}

.LogUTMS <- function(LogXDoubleFactorial, splits) {

  splits <- splits[splits > 0L]
  totalTips <- sum(splits)
  tipsMinusLengthSplits <- totalTips - length(splits)

  sum(LogXDoubleFactorial(totalTips + totalTips - 5L),
      LogXDoubleFactorial(splits + splits - 3L)) -
      LogXDoubleFactorial(tipsMinusLengthSplits + tipsMinusLengthSplits - 1L)
}

#' @rdname UnrootedTreesMatchingSplit
#' @export
LnUnrootedTreesMatchingSplit <- function(...) {
  .LogUTMS(LnDoubleFactorial, c(...))
}

#' @rdname UnrootedTreesMatchingSplit
#' @export
Log2UnrootedTreesMatchingSplit <- function(...) {
  .LogUTMS(Log2DoubleFactorial, c(...))
}
