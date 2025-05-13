#' Renumber a tree's nodes and tips
#'
#' `Renumber()` numbers the nodes and tips in a tree to conform with the
#' `phylo` standards.
#'
#' The \pkg{ape} class `phylo` is not formally defined, but expects trees' internal
#' representation to conform to certain principles: for example, nodes should
#' be numbered sequentially, with values increasing away from the root.
#'
#' `Renumber()` attempts to reformat any tree into a representation that will
#' not cause \pkg{ape} functions to produce unwanted results or to crash R.
#'
#' @template treeParam
#'
#' @examples
#' tree <- RandomTree(letters[1:10])
#' Renumber(tree)
#'
#' @return `Renumber()` returns a tree of class `phylo`, numbered in a
#' [Cladewise] fashion consistent with the expectations of \pkg{ape} functions.
#'
#' @seealso `Preorder()` provides a faster and simpler alternative, but also
#' rotates nodes.
#'
#' @template MRS
#' @family tree manipulation
#' @export
Renumber <- function(tree) {
  tree   <- ApePostorder(tree)
  edge   <- tree[["edge"]]
  nTip   <- NTip(tree)
  parent <- edge[, 1L]
  child  <- edge[, 2L]
  NODES  <- child > nTip
  TIPS   <- !NODES
  nNode  <- sum(NODES) + 1 # Root node has no edge leading to it, so add 1
  nodeLb <- tree[["node.label"]]

  tip <- child[TIPS]
  tree[["tip.label"]] <- tree[["tip.label"]][tip]
  child[TIPS] <- seq_len(nTip)

  old.node.number <- unique(parent)
  new.node.number <- rev(nTip + seq_along(old.node.number))
  renumbering.schema <- integer(nNode)
  renumbering.schema[old.node.number - nTip] <- new.node.number
  child[NODES] <- renumbering.schema[child[NODES] - nTip]
  nodeseq <- seq_len(nNode) * 2L
  parent <- renumbering.schema[parent - nTip]
  if (!is.null(nodeLb)) {
    tree[["node.label"]][new.node.number - nTip] <- 
      nodeLb[old.node.number - nTip]
  }

  tree[["edge"]][, 1] <- parent
  tree[["edge"]][, 2] <- child
  Cladewise(tree)
}

#' Generate trivial trees
#'
#' `SingleTaxonTree()` creates a phylogenetic "tree" that contains a single
#' taxon. 
#' `ZeroTaxonTree()` creates an empty `phylo` object with zero leaves or edges.
#'
#' @param label a character vector specifying the label of the tip.
#' @param lengths a numeric vector specifying the edge lengths of the tree.
#' @return `SingleTaxonTree()` returns a \code{phylo} object containing a single
#' tip with the specified label.
#'
#' @examples
#' SingleTaxonTree("Homo_sapiens")
#' plot(SingleTaxonTree("root") + BalancedTree(4))
#'
#' @keywords tree
#' @family tree manipulation
#' @family tree generation functions
#' @name TrivialTree
NULL

#' @rdname TrivialTree
#' @export
SingleTaxonTree <- function(label = "t1", lengths = NULL) {
  if (is.null(lengths)) {
    structure(list(edge = matrix(c(2L, 1L), 1, 2), tip.label = label,
                   Nnode = 1L),
              class = "phylo")
  } else {
    structure(list(edge = matrix(c(2L, 1L), 1, 2), tip.label = label,
                   Nnode = 1L, edge.length = lengths[[1]]),
              class = "phylo")
  }
}

#' @rdname TrivialTree
#' @return `ZeroTaxonTree()` returns an empty \code{phylo} object.
#' @examples
#' ZeroTaxonTree()
#' @export
ZeroTaxonTree <- function() {
  structure(list(
    # Order is consistent with ape::read.tree (but not ape::rtree...)
    edge = structure(numeric(0), dim = c(0L, 2L)),
    Nnode = 0, tip.label = character(0)
  ), class = "phylo")
}

#' Extract a subtree
#'
#' `Subtree()` safely extracts a clade from a phylogenetic tree.
#'
#' Modified from the \pkg{ape} function \code{\link[ape]{extract.clade}}, which
#' sometimes behaves unpredictably.
#' Unlike extract.clade, this function supports the extraction of "clades"
#' that constitute a single tip.
#'
#' @template preorderTreeParam
#' @param node The number of the node at the base of the clade to be extracted.
#'
#' @return `Subtree()` returns a tree of class \code{phylo} that represents a
#' clade extracted from the original tree.
#'
#' @examples
#' tree <- Preorder(BalancedTree(8))
#' plot(tree)
#' ape::nodelabels()
#' ape::nodelabels(13, 13, bg="yellow")
#'
#' plot(Subtree(tree, 13))
#'
#' @template MRS
#' @family tree manipulation
#' @export
Subtree <- function(tree, node) {
  if (is.null(treeOrder <- attr(tree, "order")) || treeOrder != "preorder") {
    stop("Tree must be in preorder")
  }
  tipLabel <- tree[["tip.label"]]
  nTip <- length(tipLabel)
  if (node <= nTip) {
    return(SingleTaxonTree(tipLabel[node]))
  }
  if (node == nTip + 1L) {
    return(tree)
  }

  edge <- tree[["edge"]]
  parent <- edge[, 1]
  child <- edge[, 2]
  subtreeParentEdge <- match(node, child)
  keepEdge <- DescendantEdges(parent, child, edge = subtreeParentEdge)
  keepEdge[subtreeParentEdge] <- FALSE

  edge <- edge[keepEdge, ]
  edge1 <- edge[, 1]
  edge2 <- edge[, 2]

  isTip <- edge2 <= nTip
  tips  <- edge2[isTip]
  new.nTip <- length(tips)
  name <- character(new.nTip)
  # method="radix" typically a few % faster than "auto"
  tipOrder <- order(tips, method = "radix")
  name[tipOrder] <- tipLabel[tips]
  edge2[isTip] <- tipOrder

  ## renumber nodes:
  nodeAdjust <- new.nTip + 1 - node
  keptLabels <- c(node, edge2[!isTip]) - nTip
  edge2[!isTip] <- edge2[!isTip] + nodeAdjust
  edge[, 1] <- edge1 + nodeAdjust
  edge[, 2] <- edge2

  ret <- structure(list(
    tip.label = name,
    Nnode = dim(edge)[[1]] - new.nTip + 1L,
    edge = edge
  ), class = "phylo", order = "preorder")
  if (!is.null(tree[["node.label"]])) {
    ret[["node.label"]] <- tree[["node.label"]][keptLabels]
  }
  
  # Return:
  ret
}

#' List ancestors
#'
#' `ListAncestors()` reports all ancestors of a given node.
#'
#' Note that if `node = NULL`, the tree's edges must be listed such that each
#' internal node (except the root) is listed as a child before it is listed
#' as a parent, i.e. its index in `child` is less than its index in `parent`.
#' This will be true of trees listed in [Preorder].
#'
#'
#' @template treeParent
#' @template treeChild
#' @param node Integer giving the index of the node or tip whose ancestors are
#'  required, or `NULL` to return ancestors of all nodes.
#'
#' @return
#' If `node = NULL`, `ListAncestors()` returns a list. Each entry _i_ contains
#' a vector containing, in order, the nodes encountered when traversing the tree
#' from node _i_ to the root node.
#' The last entry of each member of the list is therefore the root node,
#' with the exception of the entry for the root node itself, which is a
#' zero-length integer.
#'
#' If `node` is an integer, `ListAncestors()` returns a vector of the numbers of
#' the nodes ancestral to the given \code{node}, including the root node.
#' @template MRS
#'
#' @seealso
#' Implemented less efficiently in \code{phangorn:::Ancestors}, on which this
#' code is based.
#'
#' @examples
#' tree <- PectinateTree(5)
#' edge <- tree[["edge"]]
#'
#' # Identify desired node with:
#' plot(tree)
#' nodelabels()
#' tiplabels()
#'
#' # Ancestors of specific nodes:
#' ListAncestors(edge[, 1], edge[, 2], 4L)
#' ListAncestors(edge[, 1], edge[, 2], 8L)
#'
#' # Ancestors of each node, if tree numbering system is uncertain:
#' lapply(seq_len(max(edge)), ListAncestors,
#'        parent = edge[, 1], child = edge[, 2])
#'
#' # Ancestors of each node, if tree is in preorder:
#' ListAncestors(edge[, 1], edge[, 2])
#'
#' @family tree navigation
#' @export
ListAncestors <- function(parent, child, node = NULL) {
  if (is.null(node)) {
    AllAncestors(parent, child)
  } else {
    pvector <- integer(max(parent))
    pvector[child] <- parent
    anc <- function(pvector, node) {
      res <- integer(0)
      repeat {
        anc <- pvector[node]
        if (anc == 0)
          break
        res <- c(res, anc)
        node <- anc
      }
      res
    }

    # Return:
    anc(pvector, node)
  }
}

#' @describeIn ListAncestors Alias for `ListAncestors(node = NULL)`.
#'
#' @examples
#' # Alias:
#' AllAncestors(edge[, 1], edge[, 2])
#'
#' @family tree navigation
#' @export
AllAncestors <- function(parent, child) {
  res <- lapply(integer(max(parent)), integer)
  for (i in seq_along(parent)) {
    pa <- parent[i]
    res[[child[i]]] <- c(pa, res[[pa]])
  }
  res
}

#' Clade sizes
#'
#' `CladeSizes()` reports the number of nodes in each clade in a tree.
#'
#' @template treeParam
#' @param internal Logical specifying whether internal nodes should be counted
#' towards the size of each clade.
#' @param nodes Integer specifying indices of nodes at the base of clades
#' whose sizes should be returned.
#' If unspecified, counts will be provided for all nodes (including leaves).
#'
#' @return `CladeSizes()` returns the number of nodes (including leaves) that
#' are descended from each node, not including the node itself.
#'
#' @examples
#' tree <- BalancedTree(6)
#' plot(tree)
#' ape::nodelabels()
#' CladeSizes(tree, nodes = c(1, 8, 9))
#'
#' @family tree navigation
#' @export
CladeSizes <- function(tree, internal = FALSE, nodes = NULL) {
  if (length(internal) > 1 || !is.logical(internal)) {
    warning("`internal` should be a single logical value.")
    internal <- isTRUE(internal)
  }

  nTip <- NTip(tree)
  edge <- tree[["edge"]]

  size <- c(rep.int(1L, nTip), rep.int(internal, max(edge[, 1]) - nTip))
  for (i in postorder_order(tree[["edge"]])) {
    edgeI <- edge[i, ]
    parent <- edgeI[1]
    child <- edgeI[2]
    size[parent] <- size[parent] + size[child]
  }


  # Return:
  (if (is.null(nodes)) size else size[nodes]) - internal
}
