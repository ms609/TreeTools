
#' Ancestral edge
#'
#' @param edge Number of an edge
#' @template treeParent
#' @template treeChild
#' @return `AncestorEdge` returns a logical vector identifying whether each edge
#' is the immediate ancestor of the given edge.
#' @examples
#' tree <- BalancedTree(6)
#' parent <- tree$edge[, 1]
#' child <- tree$edge[, 2]
#' plot(tree)
#' ape::edgelabels()
#' AncestorEdge(5, parent, child)
#' which(AncestorEdge(5, parent, child))
#'
#' @keywords internal
#' @family tree navigation
#' @export
AncestorEdge <- function(edge, parent, child) child == parent[edge]

#' Ancestors of an edge
#'
#' Quickly identify edges that are "ancestral" to a particular edge in a tree.
#'
#' @param edge Integer specifying the number of the edge whose child edges
#' should be returned.
#' @template treeParent
#' @template treeChild
#' @param stopAt Integer or logical vector specifying the edge(s) at which to
#' terminate the search; defaults to the edges with the smallest parent,
#' which will be the root edges if nodes are numbered [Cladewise] or in
#' [Preorder].
#'
#' @return `EdgeAncestry()` returns a logical vector stating whether each edge
#' in turn is a descendant of the specified edge.
#'
#' @examples
#' tree <- PectinateTree(6)
#' plot(tree)
#' ape::edgelabels()
#' parent <- tree$edge[, 1]
#' child <- tree$edge[, 2]
#' EdgeAncestry(7, parent, child)
#' which(EdgeAncestry(7, parent, child, stopAt = 4))
#'
#' @template MRS
#' @family tree navigation
#' @export
EdgeAncestry <- function(edge, parent, child,
                         stopAt = (parent == min(parent))) {
  ret <- edge <- AncestorEdge(edge, parent, child)
  if (any(ret)) repeat {
    if (any(ret[stopAt])) return(ret)
    ret[edge <- AncestorEdge(edge, parent, child)] <- TRUE
  }
  # Return:
  ret
}

#' Most recent common ancestor
#'
#' `MRCA()` calculates the last common ancestor of specified nodes.
#'
#' `MRCA()` requires that node values within a tree increase away from the root,
#' which will be true of trees listed in `Preorder`.
#' No warnings will be given if trees do not fulfil this requirement.
#'
#' @param x1,x2 Integer specifying index of leaves or nodes whose most
#' recent common ancestor should be found.
#' @param ancestors List of ancestors for each node in a tree. Perhaps
#' produced by [`ListAncestors()`].
#'
#' @return `MRCA()` returns an integer specifying the node number of the last
#' common ancestor of `x1` and `x2`.
#'
#' @family tree navigation
#' @template MRS
#'
#' @examples
#' tree <- BalancedTree(7)
#'
#' # Verify that node numbering increases away from root
#' plot(tree)
#' nodelabels()
#'
#' # ListAncestors expects a tree in Preorder
#' tree <- Preorder(tree)
#' edge <- tree$edge
#' ancestors <- ListAncestors(edge[, 1], edge[, 2])
#' MRCA(1, 4, ancestors)
#'
#' # If a tree must be in postorder, use:
#' tree <- Postorder(tree)
#' edge <- tree$edge
#' ancestors <- lapply(seq_len(max(edge)), ListAncestors,
#'                     parent = edge[, 1], child = edge[, 2])

#'
#' @export
MRCA <- function(x1, x2, ancestors) {
  if (x1 == x2) {
    x1
  } else {
    anc1 <- ancestors[[x1]]
    anc2 <- ancestors[[x2]]
    max(intersect(anc1, anc2))
  }
}
