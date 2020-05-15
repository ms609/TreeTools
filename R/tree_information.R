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

#' Cladistic information content
#'
#' Calculate the cladistic (phylogenetic) information content of a
#' phylogenetic object, _sensu_ Thorley _et al._ (1998).
#'
#' The CIC is the logarithm of the number of binary trees that include the
#' specified topology.  A base two logarithm gives an information content in
#' bits.
#'
#' The CIC was originally proposed by Rohlf (1982), and formalised as CIC, with
#' an information-theoretic justification, by Thorley _et al_. (1998).
#' Steel and Penny (2006) term the equivalent quantity 'phylogenetic information
#' content' in the context of individual characters.
#'
#' The number of binary trees consistent with a cladogram provides a more
#' satisfactory measure of the resolution of a tree than simply
#' counting the number of edges resolved (Page, 1992).
#'
#' @param x Tree of class `phylo`, or a list thereof.
#'
#' @return Returns the phylogenetic or clustering information content
#' of the input tree(s).
#'
#' @references
#'
#' \insertRef{Page1992}{TreeTools}
#'
#' \insertRef{Rohlf1982}{TreeTools}
#'
#' \insertRef{Steel2006}{TreeTools}
#'
#' \insertRef{Thorley1998}{TreeTools}
#'
#' @family tree information functions
#' @family tree characterization functions
#'
#' @template MRS
#' @export
CladisticInfo <- function (x) UseMethod('CladisticInfo')

#' @rdname CladisticInfo
#' @export
PhylogeneticInfo <- function (x) {
  .Deprecated('CladisticInfo()')
  UseMethod('CladisticInfo')
}

CladisticInfo.phylo <- function (x) {
  Log2Unrooted(NTip(x)) - Log2TreesMatchingTree(x)
}

#' @export
CladisticInfo.list <- function (x) vapply(x, CladisticInfo, 0)
#' @export
CladisticInfo.multiPhylo <- CladisticInfo.list


#' @rdname CladisticInfo
#' @export
PhylogeneticInformation <- PhylogeneticInfo
#' @rdname CladisticInfo
#' @export
CladisticInformation <- CladisticInfo
