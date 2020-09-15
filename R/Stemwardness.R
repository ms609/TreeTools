#' 'Stemwardness' of a leaf
#'
#' Functions to describe the position of a leaf relative to the root.
#' 'Stemmier' leaves ought to exhibit a smaller root-node distance and a
#' larger sister size,
#'
#' `RootNodeDistance()` calculates the number of nodes between the chosen leaf
#' and the root of `tree`.
#' This is an unsatisfactory measure, as the range of possible
#' distances is a function of the shape of the tree.
#' As an example, leaf _X1_ in the tree `(.,(.,(.,(.,(X1,(a,b))))))`
#' falls outside the clade _(a, b)_ and has a root-node distance of 4,
#' whereas leaf _X2_ in the tree `(.,((.,(.,.)),(b,(X2,a))))`
#' falls within the clade _(a, b)_, so should be considered more 'crownwards',
#' yet has a smaller root-node distance (3).
#'
#' \insertFig{Stemwardness.png}{TreeTools}{
#'   par(mfrow = c(1, 2), mar = rep(0.3, 4))
#'   plot(ape::read.tree(text="(.,(.,(.,(.,(X1,(a,b))))));"))
#'   ape::nodelabels(1:4, 9:12)
#'   ape::edgelabels(1:2, 11:12)
#'
#'   plot(ape::read.tree(text="(.,((.,(.,.)),(b,(X2,a))));"))
#'   ape::nodelabels(1:3, c(9, 12, 13))
#'   ape::edgelabels(1, 12)
#' }
#'
#' `SisterSize()` measures the number of leaves in the clade that is sister to
#' the chosen leaf.  In the examples above, _X1_ has a sister size of 2 leaves,
#' whereas _X2_, which is 'more crownwards', has a smaller sister size (1 leaf),
#' as desired.
#'
#'
#' @template treeParam
#' @param tip Either a numeric specifying the index of a single tip, or a
#' character specifying its label.
#'
#' @return `SisterSize()` returns an integer specifying the number of leaves
#' in the clade that is sister to `tip`.
#' `RootNodeDist()` returns an integer specifying the number of nodes between
#' `tip` and the root node of `tree`.
#'
#' @references \insertRef{Asher2020}{TreeTools}
#' @examples
#' bal8 <- BalancedTree(8)
#' pec8 <- PectinateTree(8)
#'
#' SisterSize(bal8, 3)
#' SisterSize(pec8, 't3')
#' SisterSize(RootTree(pec8, 't3'), 't3')
#'
#' RootNodeDist(bal8, 3)
#' RootNodeDist(pec8, 't3')
#' RootNodeDist(RootTree(pec8, 't3'), 't3')
#' @template MRS
#' @family tree characterization functions
#' @export
SisterSize <- function (tree, tip) UseMethod('SisterSize', tip)

#' @rdname SisterSize
#' @export
SisterSize.numeric <- function (tree, tip) {
  edge <- tree$edge
  parent <- edge[edge[, 2] == tip, 1]

  # Return:
  CladeSizes(tree, nodes = parent) - 1L # Subtract the tip itself
}

#' @rdname SisterSize
#' @export
SisterSize.character <- function (tree, tip) {
  SisterSize(tree, which(tree$tip.label == tip))
}


#' @rdname SisterSize
#' @export
RootNodeDistance <- function (tree, tip) UseMethod('RootNodeDistance', tip)

#' @rdname SisterSize
#' @export
RootNodeDistance.numeric <- function (tree, tip) {
  edge <- tree$edge
  parent <- edge[edge[, 2] == tip, 1]
  depths <- NodeDepth(tree)

  # Return:
  depths[NTip(tree) + 1L] - depths[parent]

}

#' @rdname SisterSize
#' @export
RootNodeDistance.character <- function (tree, tip) {
  RootNodeDistance(tree, which(tree$tip.label == tip))
}

#' @aliases RootNodeDistance
#' @export
RootNodeDist <- RootNodeDistance
