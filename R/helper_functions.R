#' Quick sample
#' 
#' Faster than inbuilt sample because it avoids some checks
#' @param x vector to sample
#' @param len length of vector
#' @keywords internal
#' @export
SampleOne <- function (x, len = length(x)) x[sample.int(len, 1L, FALSE, NULL, FALSE)]

Assert <- function (statement) if (!statement) stop(deparse(statement), " is FALSE")

#' Edge list to edge matrix
#' @param edgeList tree edges in the format list(parent, child).
#' @return edges in the format expected by \code{tree$edge},
#'         where \code{tree} is a tree of class \code{phylo}.
#' @keywords internal
#' @export
ListToMatrix <- function (edgeList) matrix(c(edgeList[[1]], edgeList[[2]]), ncol=2)

#' Edge matrix to edge list
#' @param edge edges in the matrix format used by \code{tree$edge}, where \code{tree} is a tree of class \code{phylo}
#' @return tree edges in the format list(parent, child).
#' @export
MatrixToList <- function (edge) list(edge[, 1], edge[, 2])
