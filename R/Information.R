#' Number of trees matching a bipartition split
#'
#' Calculates the number of unrooted bifurcated trees that are consistent with
#' a bipartition split that divides taxa into groups of size `A` and `B`.
#'
#' @param A,B Integer specifying the number of taxa in each partition.
#'
#' @return `TreesMatchingSplit` returns a numeric specifying the number of trees
#' that are compatible with the given split.
#'
#' `LnTreesMatchingSplit` gives the natural logarithm of this number.
#'
#' @examples
#' TreesMatchingSplit(5, 6)
#' LnTreesMatchingSplit(5, 6)
#'
#' @template MRS
#'
#' @family split information functions
#' @export
TreesMatchingSplit <- function (A, B) {
  if (A == 0) NUnrooted(B) else
  if (B == 0) NUnrooted(A) else
  NRooted(A) * NRooted(B)
}

#' @describeIn TreesMatchingSplit Logarithm of the number of trees matching a split.
#' @export
LnTreesMatchingSplit <- function (A, B) {
  if (A == 0) LnUnrooted.int(B) else
  if (B == 0) LnUnrooted.int(A) else
  LnRooted.int(A) + LnRooted.int(B)
}

#' Character information content
#'
#' Calculates the phylogenetic information content of a given character.
#'
#' @param tokens Character vector specifying the tokens assigned to each taxon for
#' a character.  Example: `c(0, 0, 0, 1, 1, 1, '?', '-')`.
#'
#' Note that ambiguous tokens such as `(01)` are not supported, and should be
#' replaced with `?`.`
#'
#' @return `CharacterInformation` returns a numeric specifying the
#' phylogenetic information content of the character, in bits.
#'
#' @family split information functions
#' @template MRS
#' @export
CharacterInformation <- function (tokens) {
  tokenCounts <- table(tokens)
  # Our character splits our taxa into groups with the same token
  # ?s and -s are best ignored
  splits <- tokenCounts[!(names(tokenCounts) %in% c('?', '-'))]

  # Information content = -log2(probability)
  # Probability of a tree being consistent with our character is
  # n trees consistent with character / n trees with that many tips
  # NUnrootedMult(splits) / NUnrooted(splits)
  # As we are working with large numbers we can use logarithms
  lnP <- LnUnrootedMult(splits) - LnUnrooted(sum(splits))
  log2P <- lnP / log(2)
  information <- -log2P

  # Return:
  information
}

#' Phylogenetic information content of a split
#'
#' Calculate the phylogenetic information content (_sensu_
#' Steel & Penny, 2006) of a split, which reflects the probability that a
#' uniformly selected random tree will contain the split:
#' a split that is consistent with a smaller number of trees will have a higher
#' information content.
#'
#' `SplitInformation()` addresses bipartition splits, which correspond to
#' edges in an unrooted phylogeny; `MultiSplitInformation()` supports splits
#' that subdivide taxa into multiple partitions, which may correspond to
#' multi-state characters in a phylogenetic matrix.
#'
#' A simple way to characterise trees is to count the number of edges.
#' (Edges are almost, but not quite, equivalent to nodes.)
#' Counting edges (or nodes) provides a quick measaure of a tree's resolution,
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
#' [`SplitwiseInfo()`](https://ms609.github.io/TreeDist/reference/SplitwiseInfo.html),
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
#' [`ClusteringInfo()`](https://ms609.github.io/TreeDist/reference/ClusteringEntropy.html)
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
#'
#' @references
#'
#' - \insertRef{Steel2006}{TreeTools}
#'
#' @family split information functions
#'
#' @seealso
#'
#' Sum the phylogenetic information content of splits within a tree:
#'  [`TreeDist`](https://ms609.github.io/TreeDist/)`::`[`SplitwiseInfo()`](https://ms609.github.io/TreeDist/reference/SplitwiseInfo.html)
#'
#'
#' Sum the clustering information content of splits within a tree:
#'  [`TreeDist`](https://ms609.github.io/TreeDist/)`::`[`ClusteringInfo()`](https://ms609.github.io/TreeDist/reference/ClusteringEntropy.html)
#'
#'
#'
#' @template MRS
#' @export
SplitInformation <- function (A, B) {
  -(LnTreesMatchingSplit(A, B) - LnUnrooted.int(A + B)) / log(2)
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
MultiSplitInformation <- function (partitionSizes) {
  -(LnUnrootedMult(partitionSizes) - LnUnrooted.int(sum(partitionSizes))) / log(2)
}

#' Number of trees consistent with split
#'
#' Calculates the number of unrooted bifurcating trees consistent with the
#' specified multi-partition split, using the formula of Carter _et al_. (1990).
#'
#' @template splitsParam
#'
#' @return `UnrootedTreesMatchingSplit` returns an integer specifying the
#'  number of unrooted bifurcating trees consistent with the specified split.
#'
#'
#' @examples
#'  UnrootedTreesMatchingSplit(c(3, 5))
#'  UnrootedTreesMatchingSplit(c(3, 2, 1, 2))
#'
#' @references
#' See Theorem 2 in \insertRef{Carter1990}{TreeTools}
#'
#' @template MRS
#' @family split information functions
#' @export
UnrootedTreesMatchingSplit <- function (splits) {
  splits <- splits[splits > 0L]
  totalTips <- sum(splits)
  tipsMinusLengthSplits <- totalTips - length(splits)
  # use exp and log as it's just as fast, but less likely to overflow to Inf
  exp(sum(LnDoubleFactorial(totalTips + totalTips - 5L),
          LnDoubleFactorial(splits + splits - 3L)) -
        LnDoubleFactorial(tipsMinusLengthSplits + tipsMinusLengthSplits - 1L))
}
