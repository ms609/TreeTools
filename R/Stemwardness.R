#' 'Stemwardness' of a leaf
#'
#' `Stemwardness()` reports two measures of the position of a leaf relative
#' to the root.
#'
#' @template treeParam
#' @param tip Either a numeric specifying the index of a single tip, or a
#' character specifying its label.
#'
#' @return `Stemwardness()` returns a named vector of length two.
#' Entry `sisterSize` is an integer specifying the number of leaves in the clade
#' that is sister to `label`.
#' Entry `rootNodeDist` is an integer specifying the number of edges between
#' 'label' and the root node.
#'
#' @references \insertRef{Asher2020}{TreeTools}
#' @examples
#' bal8 <- BalancedTree(8)
#' pec8 <- PectinateTree(8)
#' Stemwardness(bal8, 3)
#' Stemwardness(pec8, 't3')
#' Stemwardness(RootTree(pec8, 't3'), 't3')
#' @template MRS
#' @family tree characterization functions
#' @export
Stemwardness <- function (tree, tip) {
  if (!is.numeric(tip)) tip <- which(tree$tip.label == tip)
  edge <- tree$edge
  parent <- edge[edge[, 2] == tip, 1]
  depths <- NodeDepth(tree)

  # Return:
  c(sisterSize = CladeSizes(tree, nodes = parent) - 1L, # Subtract the tip itself
    rootNodeDist = depths[NTip(tree) + 1L] - depths[parent])
}
