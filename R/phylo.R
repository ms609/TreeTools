#' Renumber a tree's nodes and tips
#'
#' `Renumber()` numbers the nodes and tips in a tree to conform with the
#' `phylo` standards.
#'
#' The 'ape' class `phylo` is not formally defined, but expects trees' internal
#' representation to conform to certain principles: for example, nodes should
#' be numbered sequentially, with values increasing away from the root.
#'
#' `Renumber()` attempts to reformat any tree into a representation that will
#' not cause 'ape' functions to produce unwanted results or to crash R.
#'
#' @template treeParam
#'
#' @examples
#' tree <- RandomTree(letters[1:10])
#' Renumber(tree)
#'
#' @return `Renumber()` returns a tree of class `phylo`, numbered in a
#' [Cladewise] fashion consistent with the expectations of 'ape' functions.
#'
#' @seealso `Preorder()` provides a faster and simpler alternative, but also
#' rotates nodes.
#'
#' @template MRS
#' @family tree manipulation
#' @export
Renumber <- function (tree) {
  tree   <- ApePostorder(tree)
  edge   <- tree$edge
  nTip   <- length(tree$tip.label)
  parent <- edge[, 1L]
  child  <- edge[, 2L]
  NODES  <- child > nTip
  TIPS   <- !NODES
  nNode  <- sum(NODES) + 1 # Root node has no edge leading to it, so add 1

  tip <- child[TIPS]
  name <- vector("character", length(tip))
  name[1:nTip] <- tree$tip.label[tip]
  tree$tip.label <- name
  child[TIPS] <- 1:nTip

  old.node.number <- unique(parent)
  new.node.number <- rev(nTip + seq_along(old.node.number))
  renumbering.schema <- integer(nNode)
  renumbering.schema[old.node.number - nTip] <- new.node.number
  child[NODES] <- renumbering.schema[child[NODES] - nTip]
  nodeseq <- (1L:nNode) * 2L
  parent <- renumbering.schema[parent - nTip]

  tree$edge[,1] <- parent
  tree$edge[,2] <- child
  Cladewise(tree)
}

#' Generate a single taxon tree
#'
#' `SingleTaxonTree()` creates a phylogenetic 'tree' that contains a single
#' taxon.
#'
#' @usage SingleTaxonTree(label)
#' @param  label a character vector specifying the label of the tip.
#' @return `SingleTaxonTree()` returns a \code{phylo} object containing a single
#' tip with the specified label.
#'
#' @examples
#' SingleTaxonTree('Homo_sapiens')
#' plot(SingleTaxonTree('root') + BalancedTree(4))
#'
#' @keywords  tree
#' @family tree manipulation
#' @family tree generation functions
#' @export
SingleTaxonTree <- function (label) {
  structure(list(edge=matrix(c(2L,1L), 1, 2), tip.label=label, Nnode=1L),
            class = 'phylo')
}

#' Extract a subtree
#'
#' `Subtree()` safely extracts a clade from a phylogenetic tree.
#'
#' Modified from the \pkg{ape} function \code{\link{extract.clade}}, which
#' sometimes behaves erratically.
#' Unlike extract.clade, this function supports the extraction of 'clades'
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
#' ape::nodelabels(13, 13, bg='yellow')
#'
#' plot(Subtree(tree, 13))
#'
#' @template MRS
#' @family tree manipulation
#' @export
Subtree <- function (tree, node) {
  if (is.null(treeOrder <- attr(tree, 'order')) || treeOrder != 'preorder') {
    stop("Tree must be in preorder")
  }
  tipLabel <- tree$tip.label
  nTip <- length(tipLabel)
  if (node <= nTip) return(SingleTaxonTree(tipLabel[node]))
  if (node == nTip + 1L) return(tree)

  edge <- tree$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  subtreeParentEdge <- match(node, child)
  keepEdge <- DescendantEdges(subtreeParentEdge, parent, child)
  keepEdge[subtreeParentEdge] <- FALSE

  edge <- edge[keepEdge, ]
  edge1 <- edge[, 1]
  edge2 <- edge[, 2]

  isTip <- edge2 <= nTip
  tips  <- edge2[isTip]
  new.nTip <- length(tips)
  name <- character(new.nTip)
  # method='radix' typically a few % faster than 'auto'
  tipOrder <- order(tips, method = 'radix')
  name[tipOrder] <- tipLabel[tips]
  edge2[isTip] <- tipOrder

  ## renumber nodes:
  nodeAdjust <- new.nTip + 1 - node
  edge2[!isTip] <- edge2[!isTip] + nodeAdjust
  edge[, 1] <- edge1 + nodeAdjust
  edge[, 2] <- edge2

  # Return:
  structure(list(
    tip.label = name,
    Nnode = dim(edge)[1] - new.nTip + 1L,
    edge = edge
  ), class = 'phylo', order = 'preorder')
}

#' Add a tip to a phylogenetic tree
#'
#' `AddTip()` adds a tip to a phylogenetic tree at a specified location.
#'
#' `AddTip()` extends \code{\link{bind.tree}}, which cannot handle
#'   single-taxon trees.
#'
#' @template treeParam
#' @param where The node or tip that should form the sister taxon to the new
#' node.  To add a new tip at the root, use `where = 0`.  By default, the
#' new tip is added to a random edge.
#' @param label Character string providing the label to apply to the new tip.
#' @param nTip,nNode,rootNode Optional integer vectors specifying number of tips and
#' nodes in `tree`, and index of root node.
#' Not checked for correctness: specifying values here trades code safety for a
#' nominal speed increase.
#'
#' @return `AddTip()` returns a tree of class \code{phylo} with an additional tip
#' at the desired location.
#'
#' @template MRS
#'
#' @seealso Add one tree to another: \code{\link{bind.tree}()}
#'
#' @examples
#' plot(tree <- BalancedTree(10))
#' ape::nodelabels()
#' ape::nodelabels(15, 15, bg='green')
#'
#' plot(AddTip(tree, 15, 'NEW_TIP'))
#'
#' @keywords tree
#' @family tree manipulation
#'
#' @export
AddTip <- function (tree,
                    where = sample.int(tree$Nnode * 2 + 2L, size = 1) - 1L,
                    label = "New tip",
                    nTip = NTip(tree),
                    nNode = tree$Nnode,
                    rootNode = RootNode(tree)
                    ) {
  if (where < 1L) where <- nTip + 1L
  newTipNumber <- nTip + 1L
  treeEdge <- tree$edge

  ## find the row of 'where' before renumbering
  if (where == rootNode) {
    case <- 1L
  } else {
      insertionEdge <- which(treeEdge[, 2] == where)
      case <- if (where <= nTip) 2L else 3L
  }
  # case = 1 -> y is bound on the root of x
  # case = 2 -> y is bound on a tip of x
  # case = 3 -> y is bound on a node of x

  # Because in all situations internal nodes need to be
  # renumbered, they are changed to negatives first, and
  # nodes eventually added will be numbered sequentially
  nodes <- treeEdge > nTip
  treeEdge[nodes] <- nTip - treeEdge[nodes]  # -1, ..., -nTip
  nextNode <- -nNode - 1L
  rootNode <- nTip - rootNode

  switch(case, { # case = 1 -> y is bound on the root of x
      treeEdge <- rbind(c(nextNode, treeEdge[1]), treeEdge, c(nextNode, newTipNumber))
      rootNode <- nextNode
    }, { # case = 2 -> y is bound on a tip of x
      beforeInsertion <- seq_len(insertionEdge)
      treeEdge[insertionEdge, 2] <- nextNode
      treeEdge <- rbind(treeEdge[beforeInsertion, ],
                        c(nextNode, where),
                        c(nextNode, newTipNumber),
                        treeEdge[-beforeInsertion, ])
    }, { # case = 3 -> y is bound on a node of x
      beforeInsertion <- seq_len(insertionEdge)

      treeEdge <- rbind(treeEdge[beforeInsertion, ],
                        c(nextNode, treeEdge[insertionEdge, 2]),
                        treeEdge[-beforeInsertion, ])
      treeEdge[insertionEdge, 2] <- nextNode

      insertionEdge <- insertionEdge + 1L
      beforeInsertion <- seq_len(insertionEdge)
      treeEdge <- rbind(treeEdge[beforeInsertion, ],
                        c(nextNode, newTipNumber),
                        treeEdge[-beforeInsertion, ])
    }
  )
  tree$tip.label <- c(tree$tip.label, label)

  nNode <- nNode + 1L
  tree$Nnode <- nNode

  ## renumber nodes:
  newNumbering <- integer(nNode)
  newNumbering[-rootNode] <- newTipNumber + 1L
  childNodes <- treeEdge[, 2] < 0L

  ## executed from right to left, so newNb is modified before x$edge:
  treeEdge[childNodes, 2] <-
    newNumbering[-treeEdge[childNodes, 2]] <-
    newTipNumber + 2:nNode
  treeEdge[, 1] <- newNumbering[-treeEdge[, 1]]

  tree$edge <- treeEdge

  # Return:
  tree
}

#' @rdname AddTip
#'
#' @details `AddTipEverywhere()` adds a tip to each edge in turn.
#'
#' @param includeRoot Logical; if `TRUE`, each position adjacent
#' to the root edge is considered to represent distinct edges; if `FALSE`,
#' they are treated as a single edge.
#' @return `AddTipEverywhere()` returns a list of class `multiPhylo` containing
#' the trees produced by adding `label` to each edge of `tree` in turn.
#'
#' @examples
#' oldPar <- par(mfrow=c(2, 4), mar=rep(0.3, 4), cex=0.9)
#'
#' backbone <- BalancedTree(4)
#' # Treating the position of the root as instructive:
#' additions <- AddTipEverywhere(backbone, includeRoot = TRUE)
#' xx <- lapply(additions, plot)
#'
#' par(mfrow=c(2, 3))
#' # Don't treat root edges as distinct:
#' additions <- AddTipEverywhere(backbone, includeRoot = FALSE)
#' xx <- lapply(additions, plot)
#'
#' par(oldPar)
#'
#' @importFrom ape is.rooted
#' @export
AddTipEverywhere <- function (tree, label = 'New tip', includeRoot = FALSE) {
  nTip <- NTip(tree)
  whichNodes <- seq_len(nTip + tree$Nnode)
  edge <- tree$edge
  root <- RootNode(edge)
  if (!includeRoot) {
    parent <- edge[, 1]
    child <- edge[, 2]
    rootChildren <- child[parent == root]

    whichNodes <- if (length(rootChildren) == 2L) {
      rootChildrenNodes <- rootChildren[rootChildren > nTip]
      whichNodes[-c(root, rootChildrenNodes[1])]
    } else {
      whichNodes[-root]
    }
  }
  lapply(whichNodes, AddTip, tree = tree, label = label, nTip = nTip,
         rootNode = root)
}

#' List all ancestral nodes
#'
#' \code{AllAncestors} lists ancestors of each parent node in a tree.
#'
#' Note that the tree's edges must be listed in an order whereby each entry in
#' \code{tr$edge[n, 1]} (with the exception of the root) occurs in
#' \code{tr$edge[i < n, 2]}.
#'
#' @template treeParent
#' @template treeChild
#'
#' @examples
#'   tr <- PectinateTree(4)
#'   plot(tr)
#'   ape::tiplabels()
#'   ape::nodelabels()
#'   edge <- tr$edge
#'   AllAncestors(edge[, 1], edge[, 2])
#'
#' @return `AllAncestors()` returns a list. Entry _i_ contains a vector containing,
#' in order, the nodes encountered when traversing the tree from node _i_ to the
#' root node.
#' The last entry of each member of the list is therefore the root node,
#' with the exception of the entry for the root node itself, which is `NULL`.
#'
#' @template MRS
#' @family tree navigation
#' @export
AllAncestors <- function (parent, child) {
  res <- vector("list", max(parent))
  for (i in seq_along(parent)) {
    pa <- parent[i]
    res[[child[i]]] <- c(pa, res[[pa]])
  }
  res
}

#' List ancestors
#'
#' Reports all ancestors of a given node
#'
#' To observe the number of a node or tip, use
#' \code{plot(tree); \link[ape]{nodelabels}(); \link[ape:nodelabels]{tiplabels}();}
#'
#' @param parent the 'parent' column of the edges property of a tree of class \code{phylo};
#' @param child the 'child' column of the edges property of a tree of class \code{phylo};
#' @param node the number of the node or tip whose ancestors are required.
#'
#' @return `ListAncestors()` returns a vector of the numbers of the nodes
#' ancestral to the given \code{node}, including the root node.
#'
#' @template MRS
#'
#' @seealso
#' - \code{phangorn:::Ancestors}, a less efficient implementation on which this
#' code is based.
#'
#' @examples
#' tree   <- ape::read.tree(text='(1, (2, (3, (4, 5))));')
#' edge   <- tree$edge
#' parent <- tree$edge[, 1]
#' child  <- tree$edge[, 2]
#' ListAncestors(parent, child, 4L)
#'
#' @family tree navigation
#' @export
ListAncestors <- function (parent, child, node) {
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

#' Clade sizes
#'
#' @template treeParam
#' @param internal Logical specifying whether internal nodes should be counted
#' towards the size of each clade.
#' @param nodes Integer specifying indices of nodes whose descendant counts
#' should be returned.
#' If unspecified, counts will be provided for all nodes (including leaves).
#'
#' @return `CladeSizes()` returns the number of nodes (including leaves) that
#' are descended from each node, not including the node itself.
#'
#' @examples
#' tree <- BalancedTree(6)
#' plot(tree)
#' ape::nodelabels()
#' CladeSizes(tree, nodes = c(8, 9))
#'
#' @family tree navigation
#' @export
CladeSizes <- function (tree, internal = FALSE, nodes = NULL) {
  if (length(internal) > 1 || !is.logical(internal)) {
    warning("`internal` should be a single logical value.")
    internal <- isTRUE(internal)
  }

  edge <- Postorder(tree$edge, renumber = FALSE)
  nTip <- NTip(tree)

  size <- c(rep(1L, nTip), rep(internal, max(edge[, 1]) - nTip))
  for (i in seq_len(nrow(edge))) {
    edgeI <- edge[i, ]
    parent <- edgeI[1]
    child <- edgeI[2]
    size[parent] <- size[parent] + size[child]
  }


  # Return:
  (if (is.null(nodes)) size else size[nodes]) - internal
}
