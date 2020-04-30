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
#' specified topology.  A base two logarithm gives an information content in bits.
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
#'
#' #TODO
#' - Detailed documentation to follow
#' - CIC.Splits to follow
#' - See also Clustering Information
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
#' Unbalanced trees contain more uneven clusters, and thus contain
#' more clustering information than balanced trees.
#'
#' This is considered a defect by Rohlf (1982) and Page (1992), on the basis
#' that tree shape does not influence the amount of _cladistic_ information
#' held within a tree.  Cladistic information sees the purpose of trees as
#' depicting the evolutionary relationships of a set of taxa.  On this view,
#' a tree is either 'right' or 'wrong'.
#'
#' Clustering information represents a different concept of the information
#' within a tree; trees serve to define groups of closely related taxa.
#' On this view, an unalanced tree makes more grouping ('nesting') statements
#' than a balanced one (Adams, 1986), corresponding to more information.
#' As such, a tree can differ from the 'true' tree in small details, yet
#' still convey much 'true' information.
#'
#'
#'
#' @inheritParams PhylogeneticInfo
#' @family tree information functions
#' @template MRS
#' @examples
#' tree <- CollapseNode(BalancedTree(10), c(12:13, 19))
#' tree <- ape::read.tree(text='((a, b), (c, d), (e, (f, (g, h))));')
#' plot(tree)
#' nodelabels()
#'
#' x <- unroot(CollapseNode(BalancedTree(8), 11))
#' plot(x)
#'
#' ClusterInfo(tree)
#' ClusterInfo(BalancedTree(8))
#' ClusterInfo(PectinateTree(8))
#'
#'
#' @references
#' \insertRef{Adams1986}{TreeTools}
#'
#' \insertRef{Page1992}{TreeTools}
#'
#' \insertRef{Rohlf1982}{TreeTools}
#'
#' @export
ClusterInfo <- function (x) UseMethod('ClusterInfo')

#' @importFrom ape keep.tip root
#' @export
ClusterInfo.phylo <- function (x) {
  nTip <- NTip(x)
  if (nTip < 3L) return (0)

  edge <- x$edge
  parent <- edge[, 1]
  child <- edge[, 2]


  rootOn1 <- RootOnNode(x, parent[child == 1L], resolveRoot = is.rooted(x))
  plot(rootOn1)
  edge <- rootOn1$edge
  parent <- edge[, 1]
  child <- edge[, 2]

  clusters <- child[parent == nTip + 1L][-1]
  clusterTips <- Descendants(rootOn1, clusters, type = 'tips')
  clusterTips[[1]] <- c(1, clusterTips[[1]])
  clusterSizes <- vapply(clusterTips, length, 1L)

  message(signif(.Entropy(clusterSizes / nTip), 3), ' -- ' ,
          signif(.Entropy(clusterSizes / sum(clusterSizes)), 3))
  .Entropy(clusterSizes / nTip) +
    sum(vapply(clusterTips[clusterSizes > 2], function (tips) {
        ClusterInfo.phylo(keep.tip(rootOn1, tips))
      }, numeric(1)))
}


.Entropy <- function (p) -sum(p * log2(p))

.H <- function (...) .Entropy(c(...) / sum(...))

ClusterInfo.Splits <- function (x) {
  lx <- as.logical(x)
  nTip <- ncol(lx)
  tips <- diag(nTip) > 0L
  rownames(tips) <- seq_len(nTip)
  allSplits <- rbind(tips, !tips, lx, !lx)
  info <- 0
  available <- diag(nTip) == 0L
  #TODO DELETE:
  dimnames(available) <- list(seq_len(nTip), seq_len(nTip))

  for (i in seq_len(ncol(lx))) {
  #for (i in 1:4) {
    availableI <- available[i, ]
    infoNeeded <- log2(2L ^ sum(availableI) - 1L)
    if (infoNeeded > 0) {
      info <- info + infoNeeded
      sisterOptions <- allSplits[, available[i, ]]
      dups <- duplicated(sisterOptions, MARGIN = 1L) &
        apply(sisterOptions, 1L, any) &
        !apply(sisterOptions, 1L, all)
      which(dups)

      ourSplits <- allSplits[dups, ]
      allSplits <- allSplits[!dups, ]

      sister <- which(ourSplits[, i])
      sisterSplits <- ourSplits[sister, ]
      stillAvailable <- ourSplits[rep(sister, sum(sisterSplits)), ]
      stillAvailable[, i] <- FALSE
      available[sisterSplits, ] <- available[sisterSplits, ] & stillAvailable

      nonSister <- 3L - sister
      otherSplits <- ourSplits[nonSister, ]
      stillAvailable <- ourSplits[rep(nonSister, sum(otherSplits)), ]
      stillAvailable[, i] <- FALSE
      available[otherSplits, ] <- available[otherSplits, ] & stillAvailable

    }

  }


}


# Ensure that tree is unrooted, or root will create an extra cluster.
#' @export
ClusterInfo.matrix <- function (x) {

  depths <- NodeDepth(x, includeTips = FALSE)
  orders <- NodeOrder(x, includeAncestor = FALSE)
  rootNode <- as.character(RootNode(x))
  orders[rootNode] <- orders[rootNode] - 1L
  sum(vapply(unique(depths), function (depth) .H(orders[depths == depth]), 0))

}
