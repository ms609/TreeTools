#' Identify taxa with long branches
#' 
#' The long branch (\acronym{LB}) score \insertCite{Struck2014}{TreeTools}
#' measures the deviation of the average pairwise patristic distance of a leaf
#' from all other leaves in a tree, relative to the average leaf-to-leaf
#' distance.
#' 
#' \insertCite{Struck2014;textual}{TreeTools} proposes the standard deviation
#' of LB scores as a measure of heterogeneity that can be compared between
#' trees; and the upper quartile of LB scores as "a representative value for
#' the taxa with the longest branches".
#' 
#' @template tree(s)Param
#' @returns `LongBranch()` returns a vector giving the long branch score for 
#' each leaf in `tree`, or a list of such vectors if `tree` is a list.
#' Results are given as raw deviations, without multiplying by 100 as proposed
#' by \insertCite{Struck2014;textual}{TreeTools}.
#' @family tree properties
#' @examples
#' tree <- BalancedTree(8, lengths = c(rep(2, 4), 5:7, rep(2, 4), rep(1, 3)))
#' lb <- LongBranch(tree)
#' tree$tip.label <- paste(tree$tip.label, signif(lb, 3), sep = ": ")
#' plot(tree, tip.col = SupportColour((1 - lb) / 2), font = 2)
#' 
#' # Standard deviation of LB scores allows comparison with other trees
#' sd(lb)
#' evenLengths <- BalancedTree(8, lengths = jitter(rep(1, 14)))
#' sd(LongBranch(evenLengths))
#' 
#' # Upper quartile identifies taxa with longest branches
#' threshold <- quantile(lb, 0.75)
#' tree$tip.label[lb > threshold]
#' @template MRS
#' @export
LongBranch <- function(tree) {
  UseMethod("LongBranch")
}

#' @importFrom ape cophenetic.phylo
#' @importFrom stats cophenetic
#' @export
LongBranch.phylo <- function(tree) {
  if (is.null(tree[["edge.length"]])) {
    warning("No edge lengths specified")
    return(NULL)
  }
  
  patristic <- cophenetic(tree)
  pairs <- colSums(patristic)
  
  # Mean pairwise distance of taxon i to all other taxa in tree
  pdi <- pairs # / (length(pairs) - 1) cancels
  # Average pairwise distance across all taxa
  pda <- sum(pairs) / length(pairs) # / (length(pairs) - 1) cancels
  
  # Return:
  (pdi / pda) - 1
}

#' @export
LongBranch.multiPhylo <- function(tree) {
  lapply(tree, LongBranch.phylo)
}

#' @export
LongBranch.list <- function(tree) {
  if (all(vapply(tree, inherits, TRUE, "phylo"))) {
    LongBranch.multiPhylo(tree)
  } else {
    stop("`tree` must be a list of 'phylo' objects")
  }
}

#' @export
LongBranch.NULL <- function(tree) NULL
