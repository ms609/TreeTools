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

#' Phylogenetic information content
#' 
#' Calculate the phylogenetic information content of a phylogenetic object.
#' 
#' 
#' #TODO
#' Splits to follow
#' 
#' @param x Tree of class `phylo`, or a list thereof.
#' 
#' @return Returns the phylogenetic or clustering information content 
#' of the input tree(s).
#' 
#' @family tree information functions
#' @seealso Clustering information content: coming soon. #TODO.
#' @template MRS
#' @export
#' @name TreeInformation

#' @rdname TreeInformation
#' @export
PhylogeneticInfo <- function (x) UseMethod(PhylogeneticInfo) 

PhylogeneticInfo.phylo <- function (x) {
  Log2Unrooted(NTip(x)) - Log2TreesMatchingTree(x)
}

#' @export
PhylogeneticInfo.list <- function (x) vapply(x, PhylogeneticInfo, 0)
#' @export
PhylogeneticInfo.multiPhylo <- PhylogeneticInfo.list


#' @rdname TreeInformation
#' @export
PhylogeneticInformation <- PhylogeneticInfo