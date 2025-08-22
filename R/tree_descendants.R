#' Identify descendant edges
#'
#' `DescendantEdges()` efficiently identifies edges that are "descended" from
#' edges in a tree.
#'
#' @template treeParent
#' @template treeChild
#' @param edge Integer specifying the number of the edge whose children are
#' required (see \code{\link[ape:nodelabels]{edgelabels}()}).
#' @param node Integer specifying the number(s) of nodes whose children are
#' required.  Specify `0` to return all nodes.  If `NULL` (the default), the 
#' `edge` parameter will be used instead.
#' @param nEdge number of edges (calculated from `length(parent)` if not
#' supplied).
#' @param includeSelf Logical specifying whether to mark `edge` as its own
#' descendant.
#' @return `DescendantEdges()` returns a logical vector stating whether each
#' edge in turn is the specified edge (if `includeSelf = TRUE`)
#' or one of its descendants.
#' @examples
#' tree <- as.phylo(0, 6)
#' plot(tree)
#' desc <- DescendantEdges(tree$edge[, 1], tree$edge[, 2], edge = 5)
#' which(desc)
#' ape::edgelabels(bg = 3 + desc)
#' @family tree navigation
#' @export
DescendantEdges <- function(parent, child, edge = NULL, node = NULL,
                            nEdge = length(parent),
                            includeSelf = TRUE) {
  postorder <- PostorderOrder(cbind(parent, child))
  
  if (is.null(node)) {
    if (is.null(edge)) {
      # Dense full matrix path
      nodeDescendants <- descendant_edges(parent, child, postorder)
      entries <- pmax(0, child - min(parent) + 1)
      ret <- matrix(FALSE, nEdge, nEdge)
      ret[entries > 0, ] <- nodeDescendants[entries, ]
      if (includeSelf) {
        diag(ret) <- TRUE
      }
      ret
    } else {
      if (length(edge) == 1L) {
        descendant_edges_single(parent, child, postorder, edge,
                                include_self = includeSelf)
      } else {
        nodeDescendants <- descendant_edges(parent, child, postorder)
        entry <- child[edge] - min(parent) + 1
        ret <- ifelse(entry > 0,
                      nodeDescendants[entry, ],
                      logical(nEdge))
        if (includeSelf) {
          ret[edge] <- TRUE
        }
        ret
      }
    }
  } else {
    nodeDescendants <- descendant_edges(parent, child, postorder)
    if (length(node) == 1 && node == 0) {
      nodeDescendants
    } else {
      nodeDescendants[node - min(parent) + 1, ]
    }
  }
}

#' Identify descendant tips
#'
#' `DescendantTips()` efficiently identifies leaves (external nodes) that are
#' "descended" from edges in a tree.
#'
#' @return `DescendantTips()` returns a logical vector stating whether each
#' leaf in turn is a descendant of the specified edge.
#' @examples
#' tips <- DescendantTips(tree$edge[, 1], tree$edge[, 2], edge = 5)
#' which(tips)
#' tiplabels(bg = 3 + tips)
#' @rdname DescendantEdges
#' @export
DescendantTips <- function(parent, child, edge = NULL,
                           node = NULL,
                           nEdge = length(parent)) {
  tips <- descendant_tips(parent, child, PostorderOrder(cbind(parent, child)))
  nTip <- dim(tips)[[2]]
  if (is.null(node)) {
    if (is.null(edge)) {
      tips[child, ]
    } else {
      tips[child[edge], ]
    }
  } else {
    if (length(node) == 1 && node == 0) {
      tips
    } else {
      tips[node, ]
    }
  }
}

#' @rdname DescendantEdges
#'
#' @return `AllDescendantEdges()` is deprecated; use `DescendantEdges()`
#' instead.
#' It returns a matrix of class logical, with row _N_ specifying whether each
#' edge is a descendant of edge _N_ (or the edge itself).
#'
#' @export
AllDescendantEdges <- function(parent, child, nEdge = length(parent)) {
  .Deprecated("DescendantEdges")
  .AllDescendantEdges(parent, child, nEdge)
}

.AllDescendantEdges <- function(parent, child, nEdge = length(parent)) {
  ret <- diag(nEdge) == 1
  blankLogical <- logical(nEdge)
  allEdges <- seq_len(nEdge)
  for (edge in rev(allEdges[child > parent[1]])) {
    ret[edge, ] <- apply(ret[parent == child[edge], ], 2, any)
    ret[edge, edge] <- TRUE
  }
  
  # Return:
  ret
}

#' Count descendants for each node in a tree
#'
#' `NDescendants()` counts the number of nodes (including leaves) directly
#' descended from each node in a tree.
#'
#' @template treeParam
#'
#' @return `NDescendants()` returns an integer listing the number of direct
#' descendants (leaves or internal nodes) for each node in a tree.
#'
#' @examples
#' tree <- CollapseNode(BalancedTree(8), 12:15)
#' NDescendants(tree)
#' plot(tree)
#' nodelabels(NDescendants(tree))
#'
#' @template MRS
#' @family tree navigation
#' @export
NDescendants <- function(tree) {
  NodeOrder(tree[["edge"]], includeAncestor = FALSE)
}
