#' Number of trees containing a tree
#'
#' Calculates the number of unrooted binary trees that are consistent with
#' a tree topology on the same leaves.
#' 
#' Remember to unroot a tree first if the position of its root is arbitrary.
#'
#' @template treeParam
#'
#' @return `TreesMatchingTree()` returns a numeric specifying the number of
#' unrooted binary trees that contain all the edges present in the input tree.
#'
#' `LnTreesMatchingTree()` gives the natural logarithm of this number.
#'
#'
#' @examples
#' partiallyResolvedTree <- CollapseNode(BalancedTree(8), 12:15)
#' TreesMatchingTree(partiallyResolvedTree)
#' LnTreesMatchingTree(partiallyResolvedTree)
#' 
#' # Number of rooted trees:
#' rootedTree <- AddTip(partiallyResolvedTree, where = 0)
#' TreesMatchingTree(partiallyResolvedTree)
#'
#' @template MRS
#'
#' @family tree information functions
#' @importFrom ape unroot
#' @export
TreesMatchingTree <- function (tree) {
  prod(NUnrooted(NodeOrder(tree)))
}

#' @rdname TreesMatchingTree
#' @export
LnTreesMatchingTree <- function (tree) {
  sum(LnUnrooted(NodeOrder(tree)))
}

#' @rdname TreesMatchingTree
#' @export
Log2TreesMatchingTree <- function (tree) {
  sum(Log2Unrooted(NodeOrder(tree)))
}

#' Tree information content
#' 
#' Calculate the phylogenetic or clustering information content of a
#' phylogenetic tree.
#' 
#' 
#' #TODO
#' Detailed explanation to follow.
#' 
#' @inheritParams TreesMatchingTree
#' @family tree information functions
#' @template MRS
#' @export
#' @name TreeInformation

#' @rdname TreeInformation
#' @export
TreePhylogeneticInfo <- function (tree) {
  Log2Unrooted(NTip(tree)) - Log2TreesMatchingTree(tree)
}

#' @rdname TreeInformation
#' @export
TreePhylogeneticInformation <- TreePhylogeneticInfo