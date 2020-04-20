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
#'
#' #TODO
#' - Splits to follow
#' - Further information to follow
#' - See also Clustering Information
#' 
#' @param x Tree of class `phylo`, or a list thereof.
#' 
#' @return Returns the phylogenetic or clustering information content 
#' of the input tree(s).
#'
#' @references
#' \insertRef{Thorley1998}{TreeTools}
#'
#' @family tree information functions
#' @template MRS
#' @export
PhylogeneticInfo <- function (x) UseMethod('PhylogeneticInfo') 

PhylogeneticInfo.phylo <- function (x) {
  Log2Unrooted(NTip(x)) - Log2TreesMatchingTree(x)
}

#' @export
PhylogeneticInfo.list <- function (x) vapply(x, PhylogeneticInfo, 0)
#' @export
PhylogeneticInfo.multiPhylo <- PhylogeneticInfo.list


#' @rdname PhylogeneticInfo
#' @export
PhylogeneticInformation <- PhylogeneticInfo

#' Clustering Information
#' 
#' 
#' #TODO 
#' 
#' update TreeDist docs to refer to this function.
#' 
#' 
#' @inheritParams PhylogeneticInfo
#' @family tree information functions
#' @template MRS
#' @examples
#' tree <- CollapseNode(BalancedTree(9), 13:15)
#' plot(tree)
#' nodelabels()
#' ClusteringInfo(tree)
#' ClusteringInfo(BalancedTree(8))
#' ClusteringInfo(PectinateTree(8))
#' 
#' 
#' @export
ClusteringInfo <- function (x) UseMethod('ClusteringInfo')

ClusteringInfo.phylo <- function (x) {
  # TODO
  # This is basically a C function written in R. 
  # Once satisfied it's working... write in C!
  
  ClusteringInfo.matrix(unroot(x)$edge)
}


.Entropy <- function (p) -sum(p * log2(p))

# Ensure that tree is unrooted, or root will create an extra cluster.
ClusteringInfo.matrix <- function (x) {
  edge <- RenumberTree(x[, 1], x[, 2])
  parent <- edge[, 1] 
  child <- edge[, 2]
  maxNode <- max(parent)
  rootNode <- parent[1]
  nTip <- rootNode - 1L
  nNode <- maxNode - nTip
  nDesc <- c(rep(1L, nTip), integer(nNode))
  
  for (i in rev(seq_len(nrow(edge)))) {
    nDesc[parent[i]] <- nDesc[parent[i]] + nDesc[child[i]]
  }
  
  orders <- NodeOrder(edge, includeAncestor = FALSE)
  
  entropy <- integer(maxNode - nTip)
  internalEdge <- edge[child > nTip, ]
  nInternal <- nrow(internalEdge)
  
  nEntities <- nTip
  for (i in seq_len(nInternal) - 1L) {
    inCluster <- orders[internalEdge[nInternal - i, 2] - nTip]
    notInCluster <- nEntities - inCluster
    entropy[i + 1L] <- .Entropy(c(inCluster, notInCluster) / nEntities)
    nEntities <- nEntities - inCluster + 1L
  }
  
  message(paste(signif(entropy, 4), collapse=', '))
  sum(entropy)
  
}