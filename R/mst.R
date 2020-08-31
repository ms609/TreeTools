#' Edges of minimum spanning tree
#'
#' Calculate or plot the ends of each edge of the minimum spanning tree of a
#' distance matrix.
#'
#' @param distances Either a matrix that can be interpreted as a distance
#' matrix, or an object of class `dist`.
#' @param plot Logical specifying whether to add the minimum spanning tree
#' to an existing plot.
#' @param x,y Numeric vectors specifying the X and Y coordinates of each
#' element in `distances`.  Necessary only if `plot = TRUE`.
#' @param \dots Additional parameters to send to `[lines()]`.
#'
#' @return `MSTLines()` returns a matrix in which each row corresponds to an
#' edge of the minimum spanning tree, and each column lists the index of the
#' entry in `distances` at which the line begins and ends.
#'
#' @seealso
#' Calculate minimum spanning tree: [`ape::mst()`].
#'
#' @references
#' \insertRef{Gower1969}{TreeTools}
#'
#' @examples
#' # Corners of an almost-regular octahedron
#' points <- matrix(c(0, 0, 2, 2, 1.1, 1,
#'                    0, 2, 0, 2, 1, 1.1,
#'                    0, 0, 0, 0, 1, -1), 6)
#' distances <- dist(points)
#' MSTLines(distances)
#' plot(points[, 1:2], ann = FALSE, asp = 1)
#' MSTLines(distances, TRUE, x = points[, 1], y = points[, 2], lwd = 2)
#' @template MRS
#' @importFrom ape mst
#' @importFrom graphics lines
#' @export
MSTEdges <- function (distances, plot = FALSE, x = NULL, y = NULL, ...) {
  umst <- ape::mst(distances)
  edges <- umst == 1L
  from <- umst
  from[edges] <- unlist(apply(edges, 1, which))
  to <- t(from)
  ends <- cbind(from[lower.tri(from) & edges], to[lower.tri(to) & edges])
  if (plot) {
    apply(ends, 1, function (edge)
      lines(x[edge], y[edge], ...))
    invisible(ends)
  } else {
    ends
  }
}
