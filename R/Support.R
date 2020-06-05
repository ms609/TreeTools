#' Frequency of splits
#'
#' `SplitFrequency()` provides a simple way to count the number of times that
#' bipartition splits, as defined by a reference tree, occur in a forest of trees.
#'
#' If multiple calculations are required, some time can be saved by using the
#' constituent functions (see examples)
#'
#'
#' @param reference A tree of class `phylo`, a `Splits` object.
#' @param forest a list of trees of class `phylo`, or a `multiPhylo` object; or a
#' `Splits` object.
#'
#' @return `SplitFrequency()` returns the number of trees in `forest` that
#' contain each split in `reference`.
#' If `reference` is a tree of class `phylo`, then the sequence will correspond
#' to the order of nodes (use `ape::nodelabels()` to view).
#' Note that the three nodes at the root of the tree correspond to a single
#' split; see the example for how these might be plotted on a tree.
#'
#' @examples
#' forest <- as.phylo(c(1, 10, 10, 100, 1000), nTip = 7)
#'
#' # Simple, but means counting each split in the forest twice:
#' tree1Freqs <- SplitFrequency(forest[[1]], forest)
#' SplitFrequency(forest[[2]], forest)
#'
#' plot(forest[[1]])
#' ape::nodelabels(tree1Freqs, node = as.integer(names(tree1Freqs)))
#'
#' @template MRS
#' @export
SplitFrequency <- function(reference, forest) {
  referenceSplits <- as.Splits(reference)
  forestSplits <- as.Splits(forest,
                            tipLabels = attr(referenceSplits, 'tip.label'))

  ret <- rowSums(vapply(forestSplits,
                        function (cf) referenceSplits %in% cf,
                        logical(length(referenceSplits))))
  names(ret) <- rownames(referenceSplits)

  # Return:
  ret

}

#' @describeIn SplitFrequency Assign a unique integer to each split
#' @param tips Integer vector specifying the tips of the tree within the chosen
#' split.
#' @template treeParam
#' @param tipIndex Character vector of tip names, in a fixed order.
#' @param powersOf2 Integer vector of same length as `tipIndex`, specifying a
#' power of 2 to be associated with each tip in turn.
#' @export
SplitNumber <- function (tips, tree, tipIndex, powersOf2) { # nocov start
  .Deprecated("SplitFrequency")
  included <- tipIndex %in% tree$tip.label[tips]
  as.character(min(c(sum(powersOf2[included]), sum(powersOf2[!included]))))
}

#' @describeIn SplitFrequency Frequency of splits in a given forest of trees
#' @export
ForestSplits <- function (forest, powersOf2) {
  .Deprecated("SplitFrequency")
  if (inherits(forest, 'phylo')) forest <- structure(list(forest), class='multiPhylo')
  tipIndex <- sort(forest[[1]]$tip.label)
  nTip <- length(tipIndex)

  # Return:
  table(vapply(forest, function (tr) {
    # +2: Don't consider root node (not a node) or first node (duplicated)
    vapply(Descendants(tr, nTip + 2L + seq_len(nTip - 3L), type='tips'),
           SplitNumber, character(1), tr, tipIndex, powersOf2)
  }, character(nTip - 3L)))
}

#' @describeIn SplitFrequency Deprecated. Listed the splits in a given tree.
#' Use as.Splits instead.
#' @export
TreeSplits <- function (tree) {
  .Deprecated("as.Splits")
} # nocov end

#' Colour for node support value
#'
#' Colour value with which to display node support.
#'
#' @param support A numeric vector of values in the range 0--1.
#' @param show1 Logical specifying whether to display values of 1.
#'              A transparent white will be returned if `FALSE`.
#' @return `SupportColour()` returns a string containing the hexadecimal code
#'  for a colour picked from a diverging scale, or `red` if a value is invalid.
#' @examples
#' SupportColour(0:4 / 4, show1 = FALSE)
#' @importFrom colorspace diverge_hcl
#' @export
SupportColour <- function (support, show1=TRUE) {
  # continuousScale <- rev(colorspace::heat_hcl(101, h=c(300, 75), c.=c(35, 95), l=c(15, 90), power=c(0.8, 1.2))) # Viridis prefered
  divergingScale <- rev(diverge_hcl(101, h=c(260, 0), c=100, l=c(50, 90), power=1.0))
  ifelse(is.na(support) | support < 0 | support > 1 | support == '', 'red',
         ifelse(support == 1 & !show1, "#ffffff00", divergingScale[(support * 100) + 1L]))
}

#' @rdname SupportColour
#' @export
SupportColor <- SupportColour
