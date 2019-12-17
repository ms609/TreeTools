#' Root a phylogenetic tree
#'
#' Roots a tree on the smallest clade containing the specified tips.
#'
#' @template treeParam
#' @template outgroupTipsParam
#'
#' @return `RootTree` returns a tree of class `phylo`, rooted on the smallest
#' clade that contains the specified tips.
#'
#' @examples
#' tree <- PectinateTree(8)
#' plot(tree)
#' ape::nodelabels()
#'
#' plot(RootTree(tree, c('t6', 't7')))
#'
#' plot(RootOnNode(tree, 12))
#' plot(RootOnNode(tree, 2))
#'
#' @seealso
#' - [ape::root]
#' - [`EnforceOutgroup`]
#'
#' @template MRS
#' @importFrom phangorn Ancestors Descendants
#' @importFrom ape root
#' @export
RootTree <- function (tree, outgroupTips) {
  tipLabels <- tree$tip.label
  if (!all(outgroupTips %in% tipLabels)) {
    stop("Outgroup tips", paste(outgroupTips, collapse=', '),
         "not found in tree's tip labels.")
  }
  if (length(outgroupTips) == 1) {
    outgroup <- outgroupTips
  } else {
    tipNos <- which(tipLabels %in% outgroupTips)
    ancestry <- unlist(Ancestors(tree, tipNos))
    ancestryTable <- table(ancestry)
    lineage <- as.integer(names(ancestryTable))
    lca <- max(lineage[ancestryTable == length(outgroupTips)])
    rootNode <- length(tipLabels) + 1L
    if (lca == rootNode) {
      lca <- lineage[lineage - c(lineage[-1], 0) != -1][1] + 1L
    }
    outgroup <- Descendants(tree, lca)[[1]]
  }

  Renumber(root(tree, outgroup, resolve.root = TRUE))
}

#' @describeIn RootTree Roots a tree on a specified internal node.
#' @param node integer specifying node (internal or tip) to set as the root.
#' @param resolveRoot logical specifying whether to resolve the root node.
#'
#' @return `RootOnNode` returns a tree of class `phylo`, rooted on the requested
#' `node` and ordered in [`Preorder`].
#'
#' @importFrom ape unroot
#' @export
RootOnNode <- function (tree, node, resolveRoot = FALSE) {
  edge <- tree$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  nodeParentEdge <- child == node
  rootEdges <- !parent %in% child
  rootNode <- parent[rootEdges][1]
  rooted <- sum(rootEdges) == 2L

  if (any(nodeParentEdge)) {
    if (rooted) {
      rootChildren <- child[rootEdges]
      spareRoot <- rootChildren == min(rootChildren) # Hit tip if present
      parent[rootEdges][!spareRoot] <- rootChildren[spareRoot]
      inverters <- EdgeAncestry(which(nodeParentEdge), parent, child) |
        rootEdges
      if (resolveRoot) {
        inverters <- inverters | nodeParentEdge
        parent[rootEdges][spareRoot] <- node
        child[rootEdges][spareRoot] <- rootNode
        child[nodeParentEdge] <- rootNode
      } else {
        if (node > rootNode) inverters <- inverters | nodeParentEdge
        deletedEdge <- which(rootEdges)[spareRoot]

        parent <- parent[-deletedEdge]
        child <- child[-deletedEdge]
        inverters <- inverters[-deletedEdge]

        parent[parent > rootNode] <- parent[parent > rootNode] - 1L
        child[child > rootNode] <- child[child > rootNode] - 1L

        tree$Nnode <- tree$Nnode - 1L

      }
    } else {
      if (resolveRoot) {
        inverters <- c(EdgeAncestry(which(nodeParentEdge), parent, child), FALSE)
        newNode <- max(parent) + 1L
        parent <- c(parent, newNode)
        child <- c(child, parent[nodeParentEdge])
        parent[nodeParentEdge] <- newNode
        tree$Nnode <- 1L + tree$Nnode
      } else {
        inverters <- EdgeAncestry(which(nodeParentEdge), parent, child)
      }
    }
    tree$edge <- RenumberTree(ifelse(inverters, child, parent),
                              ifelse(inverters, parent, child))
    attr(tree, 'order') <- 'preorder'
    tree
  } else {
    # Root position is already correct
    if (rooted && !resolveRoot) {
      unroot(tree)
    } else if (!rooted && resolveRoot) {
      RootOnNode(tree, max(child[rootEdges]), resolveRoot = TRUE)
    } else {
      Preorder(tree)
    }
  }
}

#' Collapse nodes on a phylogenetic tree
#'
#' Collapses specified nodes or edges on a phylogenetic tree, resulting in
#' polytomies.
#'
#' @template treeParam
#' @param nodes,edges Integer vector specifying the nodes or edges in the tree
#'  to be dropped.
#' (Use \code{\link[ape]{nodelabels}} or \code{\link[ape:nodelabels]{edgelabels}}
#' to view numbers on a plotted tree.)
#'
#' @return `CollapseNode` and `CollapseEdge` return a tree of class `phylo`,
#' corresponding to `tree` with the specified nodes or edges collapsed.
#' The length of each dropped edge will (naively) be added to each descendant
#' edge.
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' oldPar <- par(mfrow=c(2, 1), mar=rep(0.5, 4))
#'
#' tree <- rtree(7)
#' plot(tree)
#' nodelabels()
#' edgelabels(round(tree$edge.length, 2), cex=0.6, frame='n', adj=c(1, -1))
#'
#' newTree <- CollapseNode(tree, c(12, 13))
#' plot(newTree)
#' nodelabels()
#' edgelabels(round(newTree$edge.length, 2), cex=0.6, frame='n', adj=c(1, -1))
#'
#' par(oldPar)
#'
#' @author  Martin R. Smith
#' @export
CollapseNode <- function (tree, nodes) {
  if (length(nodes) == 0) return (tree)

  edge <- tree$edge
  lengths <- tree$edge.length
  hasLengths <- !is.null(lengths)
  parent <- edge[, 1]
  child <- edge[, 2]
  root <- min(parent)
  nTip <- root - 1L
  maxNode <- max(parent)
  edgeBelow <- order(child)
  edgeBelow <- c(edgeBelow[1:(root-1L)], NA, edgeBelow[-(1:root-1L)])
  nodes <- unique(nodes)

  if (!inherits(tree, 'phylo')) stop ("tree must be an object of class phylo")
  if (!all(nodes %in% (root + 1L):maxNode)) stop("nodes must be integers between ",
                                                 root + 1L, " and ", maxNode)

  keptEdges <- -edgeBelow[nodes]

  for (node in rev(sort(nodes))) {
    newParent <- parent[edgeBelow[node]]
    if (hasLengths) lengths[parent == node] <- lengths[parent == node] + lengths[child == node]
    parent[parent == node] <- newParent
  }

  newNumber <- c(seq_len(nTip), nTip + cumsum(root:maxNode %in% parent))

  tree$edge <-cbind(newNumber[parent[keptEdges]], newNumber[child[keptEdges]])
  tree$edge.length <- lengths[keptEdges]
  tree$Nnode <- tree$Nnode - length(nodes)

  # TODO renumber nodes sequentially
  # TODO should probably re-write this in C++.
  tree
}

#' @rdname CollapseNode
#' @export
CollapseEdge <- function (tree, edges) {
  CollapseNode(tree, tree$edge[edges, 2])
}
