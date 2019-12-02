#' Reorder edges of a phylogenetic tree
#'
#' Wrappers for the C functions called by
#' \code{ape::\link[ape:reorder.phylo]{reorder.phylo}}.
#' These call the C functions directly, so are faster -- but don't perform
#' as many checks on user input.  Bad input could crash R.
#'
#'
#' @param nTip,nNode,nEdge Integer specifying the number of tips, nodes
#' and edges in the input tree.
#' @template treeParent
#' @template treeChild
#'
#' @return `NeworderPruningwise` returns an integer vector specifying the
#' pruningwise order of edges within a tree.
#'
#' @examples
#' nTip <- 8L
#' tree <- BalancedTree(nTip)
#' edge <- tree$edge
#' pruningwise <- NeworderPruningwise(nTip, tree$nNode, edge[, 1], edge[, 2],
#'                              dim(edge)[1])
#' cladewise <- NeworderPhylo(nTip,edge[, 1], edge[, 2], dim(edge)[1], 1L)
#' postorder <- NeworderPhylo(nTip,edge[, 1], edge[, 2], dim(edge)[1], 2L)
#'
#' tree$edge <- tree$edge[pruningwise, ]
#'
#' @author
#'  - C algorithm: Emmanuel Paradis
#'  - R wrapper: Martin R. Smith
#' @family C wrappers
#' @keywords internal
#' @useDynLib TreeTools, .registration = TRUE
#' @name Neworder
#' @export
NeworderPruningwise <- function (nTip, nNode, parent, child, nEdge) {
  .C('ape_neworder_pruningwise', as.integer(nTip), as.integer(nNode),
     as.integer(parent), as.integer(child), as.integer(nEdge),
     integer(nEdge))[[6]]
}

#' @rdname Neworder
#' @param `whichwise` Integer specifying whether to order edges (1)
#' cladewise; or (2) in postorder.
#' @return `NeworderPhylo` returns an integer vector specifying the order
#' of edges under the ordering sequence specified by `whichwise`.
#' @keywords internal
#' @export
NeworderPhylo <- function (nTip, parent, child, nb.edge, whichwise) {
  .C('ape_neworder_phylo', as.integer(nTip), as.integer(parent),
     as.integer(child), as.integer(nb.edge), integer(nb.edge),
     as.integer(whichwise), NAOK = TRUE)[[5]]
}

#' Renumber a tree in preorder
#'
#' Wrapper for the C function `preorder_edges_and_nodes`, which
#' renumbers internal nodes and orders edges in preorder, in an order
#' guaranteed to be identical for any tree of an equivalent topology.
#' At each node, child edges are arranged from left to right according to the
#' lowest-numbered tip in the subtree subtended by each edge.
#'
#' @template treeParent
#' @template treeChild
#'
#' @return `RenumberTree` returns an edge matrix for a tree of class `phylo`
#' following the usual preorder convention for edge and node numbering.
#' @family C wrappers
#' @keywords internal
#' @export
RenumberTree <- function (parent, child) {
  .Call(`_TreeTools_preorder_edges_and_nodes`, parent, child)
}

#' @describeIn RenumberTree Returns in list format
#' @return `RenumberEdges` returns a list whose two entries correspond
#' to the new parent and child vectors.
#' @family C wrappers
#' @keywords internal
#' @export
RenumberEdges <- function (parent, child, nEdge = length(parent)) {
  oenn <- .Call(`_TreeTools_preorder_edges_and_nodes`, parent, child)

  # Return:
  list(oenn[, 1], oenn[, 2])
}

#' Reorder trees
#'
#' A wrapper for \code{ape:::.reorder_ape}.
#' Calling this C function directly is approximately twice as fast as using
#' \code{ape::\link[ape:reorder.phylo]{cladewise}} or
#' \code{ape::\link[ape:reorder.phylo]{postorder}}
#'
#' For `Cladewise`, `Postorder` and `Pruningwise`, all nodes must be binary;
#'  [ape::collapse.singles] and [ape::multi2di] may help.
#'
#' `Preorder` is more robust: it supports polytomies, nodes can be numbered
#' in any sequence, and edges can be listed in any order in the input tree.
#' It has a unique output for each tree topology, allowing unique trees
#' to be detected by comparing sorted edge matrices alone.
#'
#' A tree in preorder is numbered starting from the root node.
#' Each node is numbered in the sequence in which it is encountered, and
#' each edge is listed in the sequence in which it is visited.
#'
#' Child edges of a node are sorted from left to right in order of the
#' smallest descendant tip; i.e. an edge leading to tip 1 will be to the left
#' of an edge leading to a subtree that contains tip 2.
#'
#' Numbering begins by following the leftmost edge of the root node,
#' and sorting its descendant subtree into preorder.
#' Then, the next edge at the root node is followed, and its descendants
#' sorted into preorder, until each edge has been visited.
#'
#' @template treeParam
#' @template nTipParam
#' @param edge (optional) the value of tree$edge
#'
#' @return A tree with nodes following the specified numbering scheme
#' @author
#'  `Preorder`: Martin R. Smith.
#'
#' `Cladewise`, `Postorder` and `Pruningwise`: modified by Martin R. Smith from
#' \code{.reorder_ape} in \pkg{ape} (Emmanuel Paradis)
#'
#' @family C wrappers
#' @keywords internal
#' @export
Cladewise <- function (tree, nTip = NULL, edge = tree$edge) {
  if (!is.null(attr(tree, "order")) && attr(tree, "order") == "cladewise") {
    return(tree)
  }
  if (is.null(nTip)) nTip <- length(tree$tip.label)
  if (is.null(edge)) edge <- tree$edge
  nb.edge <- dim(edge)[1]
  nb.node <- tree$Nnode
  if (nb.node == 1) return(tree)
  if (nb.node >= nTip) stop("`tree` apparently badly conformed")

  neworder <- NeworderPhylo(nTip, edge[, 1], edge[, 2], nb.edge, 1)

  tree$edge <- edge[neworder, ]
  if (!is.null(tree$edge.length)) tree$edge.length <- tree$edge.length[neworder]
  attr(tree, "order") <- "cladewise"
  tree
}

#' @describeIn Cladewise Reorder tree in Postorder
#' @export
Postorder <- function (tree, nTip = length(tree$tip.label), edge = tree$edge) {
  if (!is.null(attr(tree, "order")) && attr(tree, "order") == "postorder") {
    return(tree)
  }
  nb.edge <- dim(edge)[1]
  nb.node <- tree$Nnode
  if (nb.node == 1) return(tree)
  if (nb.node >= nTip) stop("`tree` apparently badly conformed")
  neworder <- NeworderPhylo(nTip, edge[, 1], edge[, 2], nb.edge, 2)
  tree$edge <- edge[neworder, ]
  if (!is.null(tree$edge.length)) tree$edge.length <- tree$edge.length[neworder]
  attr(tree, "order") <- "postorder"
  tree
}

#' @describeIn Cladewise Reorder parent and child edges in Postorder
#' @template treeParent
#' @template treeChild
#' @template nTipParam
#' @export
PostorderEdges <- function (parent, child,
                            nEdge = length(parent),
                            nNode = nEdge / 2L,
                            nTip = nNode + 1L) {
  if (nNode == 1) return(list(parent, child))
  newOrder <- NeworderPhylo(nTip, parent, child, nEdge, 2)
  list(parent[newOrder], child[newOrder])
}

#' @describeIn Cladewise Reorder tree Pruningwise
#' @export
Pruningwise <- function (tree, nTip = length(tree$tip.label),
                         edge = tree$edge) {
  if (!is.null(attr(tree, "order")) && attr(tree, "order") == 'pruningwise') {
    return(tree)
  }
  nb.edge <- dim(edge)[1]
  nb.node <- tree$Nnode
  if (nb.node == 1) return(tree)
  if (nb.node >= nTip) stop("`tree` apparently badly conformed")
  tree <- Cladewise(tree, nTip, edge)
  neworder <- NeworderPruningwise(nTip, nb.node, tree$edge[, 1],
                                  tree$edge[, 2], nb.edge)
  tree$edge <- tree$edge[neworder, ]
  if (!is.null(tree$edge.length)) tree$edge.length <- tree$edge.length[neworder]
  attr(tree, "order") <- 'pruningwise'
  tree
}

#' @describeIn Cladewise Reorder tree in Preorder (special case of cladewise)
#' @export
Preorder <- function (tree) {
  startOrder <- attr(tree, 'order')
  if (startOrder == 'preorder') {
    tree
  } else {
    edge <- tree$edge
    parent <- edge[, 1]
    child <- edge[, 2]
    tree$edge <- RenumberTree(parent, child)
    attr(tree, 'order') <- 'preorder'
    # Return:
    tree
  }
}


#' Reorder tips
#'
#' \code{RenumberTips(tree, tipOrder)} sorts the tips of a phylogenetic tree
#' such that the indices in \code{tree$edge[, 2]} correspond to the order of
#' tips given in \code{tipOrder}
#'
#' @template treeParam
#' @param tipOrder A character vector containing the values of
#'        \code{tree$tip.label} in the desired sort order
#'
#' @examples
#' data(Lobo) # Loads the phyDat object Lobo.phy
#' tree <- RandomTree(Lobo.phy)
#' tree <- RenumberTips(tree, names(Lobo.phy))
#'
#' @template MRS
#' @export
RenumberTips <- function (tree, tipOrder) {
  startOrder <- tree$tip.label
  if (identical(startOrder, tipOrder)) return (tree)
  if (length(startOrder) != length(tipOrder)) {
    stop("Tree labels and tipOrder must match")
  }

  nTip <- length(startOrder)
  child <- tree$edge[, 2]
  tips <- child <= nTip

  matchOrder <- match(startOrder, tipOrder)
  if (any(is.na(matchOrder))) stop("All tree labels must occur in tipOrder")
  tree$edge[tips, 2] <- matchOrder[tree$edge[tips, 2]]
  tree$tip.label <- tipOrder
  tree
}
