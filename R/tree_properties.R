#' Descendant Edges
#'
#' Quickly identifies edges that are 'descended' from a particular edge in a tree
#'
#' @param edge number of the edge whose child edges are required
#' @template treeParent
#' @template treeChild
#' @param nEdge number of edges (calculated from length(parent) if not supplied)
#' @return `DescendantEdges` returns a logical vector stating whether each edge in turn is a descendant of the specified edge
#'         (or the edge itself)
#' @family tree navigation
#' @export
DescendantEdges <- function (edge, parent, child, nEdge = length(parent)) {
  ret <- logical(nEdge)
  edgeSister <- match(parent[edge], parent[-edge])
  if (edgeSister >= edge) {
    # edgeSister is really 1 higher than you think, because we knocked out edge 'edge' in the match
    ret[edge:edgeSister] <- TRUE
    
    # Return:
    ret
  } else {
    nextEdge <- edge
    revParent <- rev(parent)
    repeat {
      if (revDescendant <- match(child[nextEdge], revParent, nomatch=FALSE)) {
        nextEdge <- 1 + nEdge - revDescendant
      } else break;
    }
    ret[edge:nextEdge] <- TRUE 
    
    # Return:
    ret
  }
}

#' All Descendant Edges
#'
#' @return `AllDescendantEdges` returns a matrix of class logical, with row N specifying whether each edge is a descendant of edge N
#'         (or the edge itself)
#' @describeIn DescendantEdges Quickly identifies edges that are 'descended' from each edge in a tree
#' @family tree navigation
#' @export
AllDescendantEdges <- function (parent, child, nEdge = length(parent)) {
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

#' Ancestral edge
#'
#' @param edge Number of an edge
#' @template treeParent
#' @template treeChild
#' @return a logical vector identifying whether each edge is the edge that is ancestral to the given edge.
#' @keywords internal
#' @family tree navigation
#' @export
AncestorEdge <- function (edge, parent, child) child == parent[edge]

#' EdgeAncestry
#'
#' Quickly identifies edges that are 'ancestral' to a particular edge in a tree
#'
#' @param edge number of the edge whose child edges are required
#' @template treeParent
#' @template treeChild
#' @param stopAt number of the edge at which the search should terminate; defaults to the root edges
#' @return a logical vector stating whether each edge in turn is a descendant of the specified edge
#'
#' @author Martin R. Smith
#' @family tree navigation
#' @export
EdgeAncestry <- function (edge, parent, child, stopAt = (parent==min(parent))) {
  ret <- edge <- AncestorEdge(edge, parent, child)
  repeat {
    if (any(ret[stopAt])) return(ret)
    ret[edge <- AncestorEdge(edge, parent, child)] <- TRUE    
  }
}

#' Most Recent Common Ancestor
#' @param tip1,tip2 Integer specifying index of tips whose most recent common 
#' ancestor should be found.
#' @param ancestors Output of [`AllAncestors`] for the tree in question
#' 
#' @family tree navigation
#' @author Martin R. Smith
#' @export
MRCA <- function(tip1, tip2, ancestors) {
  anc1 <- ancestors[[tip1]]
  anc2 <- ancestors[[tip2]]
  max(intersect(anc1, anc2))
}

#' Non-duplicate root
#' 
#' Identify, for each edge, whether it is not a duplicate of the root edge
#' 
#' @template treeParent
#' @template treeChild
#' @template treeNEdgeOptional
#' 
#' @author Martin R. Smith
#' @export
#' @family tree navigation
#'  
#' @keywords internal
NonDuplicateRoot <- function (parent, child, nEdge = length(parent)) {
  notDuplicateRoot <- !logical(nEdge)
  rightSide <- DescendantEdges(1, parent, child, nEdge)
  nEdgeRight <- sum(rightSide)
  if (nEdgeRight == 1) {
    notDuplicateRoot[2] <- FALSE
  } else if (nEdgeRight == 3) {
    notDuplicateRoot[4] <- FALSE
  } else {
    notDuplicateRoot[1] <- FALSE
  }
  notDuplicateRoot
}

#' Number of partitions in a tree
#' 
#' @param tree A phylogenetic tree of class `phylo`, or a list of such trees
#' (of class `list` or `multiPhylo`), or a vector of integers.
#' 
#' @return Integer specifying the number of partitions in the specified trees,
#' or in a rooted tree with `n` tips.
#' 
#' @examples {
#' NPartitions(8)
#' NPartitions(ape::rtree(8))
#' }
#' @author Martin R. Smith
#' @importFrom ape collapse.singles
#' @export
NPartitions <- function (tree) {
  if (class(tree) == 'phylo') {
    collapse.singles(tree)$Nnode - 1L - TreeIsRooted(tree)
  } else if (mode(tree) == 'list') {
    vapply(tree, NPartitions, 1L)
  } else if (mode(tree) == 'numeric') {
    tree - 3L
  } else {
    stop("tree in unsupported format.")
  }
}
#' @rdname NPartitions
#' @export
NSplits <- NPartitions

#' Is tree rooted?
#' 
#' Faster alternative to of `ape::is.rooted`.
#' 
#' @param tree A phylogenetic tree of class phylo.
#' @return Logical specifying whether a root node is resolved.
#' 
#' @author Martin R. Smith
#' @export
TreeIsRooted <- function (tree) {
  edge <- tree$edge
  parent <- edge[, 1]
  sum(parent == min(parent)) < 3L
}
