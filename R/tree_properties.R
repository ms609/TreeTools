#' Identify descendant edges
#'
#' Quickly identify edges that are 'descended' from a particular edge in a tree.
#'
#' @param edge Integer specifying the number of the edge whose child edges are
#' required (see \code{\link[ape:nodelabels]{edgelabels}}).
#' @template treeParent
#' @template treeChild
#' @param nEdge number of edges (calculated from `length(parent)` if not
#' supplied).
#' @return `DescendantEdges` returns a logical vector stating whether each edge
#' in turn is a descendant of the specified edge (or the edge itself).
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
#' @return `AllDescendantEdges` returns a matrix of class logical, with row N
#' specifying whether each edge is a descendant of edge N (or the edge itself).
#' @describeIn DescendantEdges Quickly identifies edges that are 'descended'
#' from each edge in a tree.
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
AncestorEdge <- function (edge, parent, child) child == parent[edge]

#' Ancestros of an edge
#'
#' Quickly identify edges that are 'ancestral' to a particular edge in a tree.
#'
#' @param edge Integer specifying the number of the edge whose child edges
#' should be returned.
#' @template treeParent
#' @template treeChild
#' @param stopAt Integer or logical vector specifying the edge(s) at which to
#' terminate the search; defaults to the edges with the smallest parent,
#' which will be the root edges if nodes are numbered cladewise or in Preorder.
#' @return `EdgeAncestry` returns a logical vector stating whether each edge in
#' turn is a descendant of the specified edge.
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
EdgeAncestry <- function (edge, parent, child, stopAt = (parent==min(parent))) {
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
#' What is the last common ancestor of the specified tips?
#'
#' @param tip1,tip2 Integer specifying index of tips whose most recent common
#' ancestor should be found.
#' @param ancestors Output of [`AllAncestors`] for the tree in question
#'
#' @return `MRCA` returns an integer specifying the node number of the last
#' common ancestor of `tip1` and `tip2`.
#'
#' @family tree navigation
#' @template MRS
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
#' @return A symmetrical matrix listing the number of edges that must be
#' traversed to travel from each numbered edge to each other.
#' The two edges straddling the root of a rooted tree
#' are counted as a single edge.  Add a 'root' tip using [`AddTip`] if the
#' position of the root is significant.
#'
#' @examples
#'
#' tree <- BalancedTree(5)
#' plot(tree)
#' ape::edgelabels()
#'
#' EdgeDistances(tree)
#'
#' @family tree navigation
#' @template MRS
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

  origOrder <- match(tree$edge[, 2], edge[[2]])

  # Return:
  ret[origOrder, origOrder]
}

#' Non-duplicate root
#'
#' Identify, for each edge, whether it denotes a different partition from
#' the root edge.
#' The first edge of the input tree must be a root edge; this can be
#' accomplished using `Preorder`.
#'
#' @template treeParent
#' @template treeChild
#' @template treeNEdgeOptional
#'
#' @return `NonDuplicateRoot` returns a logical vector of length `nEdge`,
#' specifying `TRUE` unless an edge identifies the same partition as
#' the root edge.
#'
#' @examples
#' tree <- Preorder(BalancedTree(8))
#' edge <- tree$edge
#' parent <- edge[, 1]
#' child <- edge[, 2]
#'
#' which(!NonDuplicateRoot(parent, child))
#'
#' @keywords internal
#' @template MRS
#' @family tree navigation
#' @export
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

#' Number of distinct partitions
#'
#' How many unique bipartition splits occur in a tree or object?
#'
#' @param x A phylogenetic tree of class `phylo`, or a list of such trees
#' (of class `list` or `multiPhylo`), or a `Splits` object,
#' or a vector of integers.
#'
#' @return `NSplits` returns an integer specifying the number of partitions in
#'  the specified objects, or in a rooted tree with `n` tips.
#'
#' @examples
#' NSplits(8L)
#' NSplits(PectinateTree(8))
#' NSplits(as.Splits(BalancedTree(8)))
#'
#'
#' @template MRS
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
NSplits.multiPhylo <- function (x) vapply(x, NSplits, numeric(1L))

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
#' @param tree A phylogenetic tree of class `phylo`.
#' @return `TreeIsRooted()` returns a logical specifying whether a root node is
#' resolved.
#'
#' @examples
#' TreeIsRooted(BalancedTree(6))
#'
#' @template MRS
#' @export
TreeIsRooted <- function (tree) {
  edge <- tree$edge
  parent <- edge[, 1]
  sum(parent == min(parent)) < 3L
}
