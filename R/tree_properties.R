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

#' Distance between edges
#'
#' Number of nodes that must be traversed to navigate from each edge to
#' each other edge within a tree
#'
#' @template treeParam
#'
#' @family tree navigation
#' @author Martin R. Smith
#' @export
EdgeDistances <- function (tree) {
  edge <- tree$edge
  nEdge <- dim(edge)[1]
  edge <- RenumberEdges(edge[, 1], edge[, 2], nEdge = nEdge)
  parent <- edge[[1]]
  child <- edge[[2]]
  ancs <- AllAncestors(parent, child)
  ret <- matrix(0L, ncol = nEdge, nrow = nEdge)
  for (i in seq_len(nEdge - 1L)) {
    ancI <- ancs[[child[i]]]
    nAncI <- length(ancI)
    for (j in seq(from = i + 1L, to = nEdge)) {
      ancJ <- ancs[[child[j]]]
      intersection <- intersect(ancI, ancJ)
      if (length(intersection) > 1L) {
        # On same side of root
        mrca <- max(intersection)
        if (child[i] %in% ancJ || child[j] %in% ancI) {
          addOne <- 0L
        } else {
          addOne <- 1L
        }
      } else {
        # On opposite sides of root
        mrca <- intersection
        addOne <- 0L
      }
      u <- union(ancI, ancJ)
      ret[i, j] <- sum(u > mrca, addOne)
    }
  }
  ret[lower.tri(ret)] <- t(ret)[lower.tri(ret)]

  # Return:
  ret
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

#' Number of distinct partitions in a tree
#'
#' @param x A phylogenetic tree of class `phylo`, or a list of such trees
#' (of class `list` or `multiPhylo`), or a `Splits` object,
#' or a vector of integers.
#'
#' @return Integer specifying the number of partitions in the specified trees,
#' or in a rooted tree with `n` tips.
#'
#' @examples {
#'   NSplits(8L)
#'   NSplits(PectinateTree(8))
#'   NSplits(as.Splits(BalancedTree(8)))
#' }
#'
#' @author Martin R. Smith
#'
#' @family Splits operations
#' @importFrom ape collapse.singles
#' @export
NSplits <- function (x) UseMethod('NSplits')

#' @rdname NSplits
#' @export
NPartitions <- NSplits

#' @rdname NSplits
#' @export
NSplits.phylo <- function (x) collapse.singles(x)$Nnode - 1L - TreeIsRooted(x)

#' @rdname NSplits
#' @export
NSplits.multiPhylo <- function (x) vapply(x, NSplits, integer(1))

#' @rdname NSplits
#' @export
NSplits.list <- NSplits.multiPhylo

#' @rdname NSplits
#' @export
NSplits.Splits <- function (x) nrow(x)

#' @rdname NSplits
#' @export
NSplits.numeric <- function (x) x - 3L


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
