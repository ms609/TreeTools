#' neworder_phylo
#'
#' Wrapper for the ape function
#' @family C wrappers
#' @keywords internal
#' @useDynLib TreeTools, .registration = TRUE
#' @export
NeworderPhylo <- function (nTip, parent, child, nb.edge, whichwise) {
  .C('ape_neworder_phylo', as.integer(nTip), as.integer(parent),
     as.integer(child), as.integer(nb.edge), integer(nb.edge),
     as.integer(whichwise), NAOK = TRUE)[[5]]
}

#' neworder_pruningwise
#'
#' Wrapper for the ape function
#' @family C wrappers
#' @keywords internal
#' @export
NeworderPruningwise <- function (nTip, nb.node, parent, child, nb.edge) {
  .C('ape_neworder_pruningwise', as.integer(nTip), as.integer(nb.node),
     as.integer(parent), as.integer(child), as.integer(nb.edge),
     integer(nb.edge))[[6]]
}

#' Renumber a tree
#'
#' Order edges and number nodes
#'
#' Wrapper for the C function `preorder_edges_and_nodes`
#'
#' @return an edge matrix for a tree in following the usual preorder convention
#'  for edge and node numbering
#' @family C wrappers
#' @keywords internal
#' @export
RenumberTree <- function (parent, child) {
  .Call(`_TreeTools_preorder_edges_and_nodes`, parent, child)
}

#' @describeIn RenumberTree Instead returns a list containing two items
#'  corresponding to the new parent and child vectors
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
#'        \code{tree$tip.label} in the desired sort order, or an object
#'        (perhaps of class `phylo` or `Splits`) with tip labels.
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
  newOrder <- .TipLabels(tipOrder)
  if (identical(startOrder, newOrder)) return (tree)
  if (length(startOrder) != length(newOrder)) {
    stop("Tree labels and tipOrder must match")
  }

  nTip <- length(startOrder)
  child <- tree$edge[, 2]
  tips <- child <= nTip

  matchOrder <- match(startOrder, newOrder)
  if (any(is.na(matchOrder))) stop("All tree labels must occur in tipOrder")
  tree$edge[tips, 2] <- matchOrder[tree$edge[tips, 2]]
  tree$tip.label <- newOrder
  tree
}
