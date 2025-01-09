#' Minimum spanning tree
#'
#' Calculate or plot the minimum spanning tree \insertCite{Gower1969}{TreeTools}
#' of a distance matrix.
#'
#' @param distances Either a matrix that can be interpreted as a distance
#' matrix, or an object of class `dist`.
#' @param plot Logical specifying whether to add the minimum spanning tree
#' to an existing plot.
#' @param x,y Numeric vectors specifying the X and Y coordinates of each
#' element in `distances`.  Necessary only if `plot = TRUE`.
#' @param \dots Additional parameters to send to `[lines()]`.
#'
#' @return `MSTEdges()` returns a matrix in which each row corresponds to an
#' edge of the minimum spanning tree, listed in non-decreasing order of length.
#' The two columns contain the indices of the entries in `distances` that
#' each edge connects, with the lower value listed first.
#'
#' @seealso
#' Slow implementation returning the association matrix of the minimum spanning
#' tree: [`ape::mst()`].
#'
#' @references \insertAllCited{}
#'
#' @examples
#' # Corners of an almost-regular octahedron
#' points <- matrix(c(0, 0, 2, 2, 1.1, 1,
#'                    0, 2, 0, 2, 1, 1.1,
#'                    0, 0, 0, 0, 1, -1), 6)
#' distances <- dist(points)
#' mst <- MSTEdges(distances)
#' MSTLength(distances, mst)
#' plot(points[, 1:2], ann = FALSE, asp = 1)
#' MSTEdges(distances, TRUE, x = points[, 1], y = points[, 2], lwd = 2)
#' @template MRS
#' @family utility functions
#' @importFrom graphics lines
#' @export
MSTEdges <- function(distances, plot = FALSE, x = NULL, y = NULL, ...) {
  ends <- MinimumSpanningTree(distances)
  if (plot) {
    apply(ends, 1, function(edge)
      lines(x[edge], y[edge], ...))
    invisible(ends)
  } else {
    ends
  }
}

#' @rdname MSTEdges
#' @param mst Optional parameter specifying the minimum spanning tree in the
#' format returned by `MSTEdges()`; if `NULL`, calculated from `distances`.
#' @return `MSTLength()` returns the length of the minimum spanning tree.
#' @export
MSTLength <- function(distances, mst = NULL) {
  distMat <- as.matrix(distances)
  if (is.null(mst)) mst <- MSTEdges(distances)
  sum(apply(mst, 1L, function(x) distMat[x[1], x[2]]))
}

MinimumSpanningTree <- function(distances) UseMethod("MinimumSpanningTree")

#' @export
MinimumSpanningTree.dist <- function(distances) {
  minimum_spanning_tree(order(distances, decreasing = TRUE) - 1L)
}

#' @export
MinimumSpanningTree.matrix <- function(distances) {
  dims <- dim(distances)
  if (dims[1] != dims[2]) {
    stop("`distances` must be a square matrix.")
  }
  dists <- distances[lower.tri(distances)]
  minimum_spanning_tree(order(dists, decreasing = TRUE) - 1L)
}
