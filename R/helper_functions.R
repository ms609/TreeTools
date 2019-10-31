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

#' Add Tree to Start of List
#'
#' `UnshiftTree` adds a phylogenetic tree to the start of a list of trees.
#' This is useful where the class of a list of trees is unknown.
#' Adding a tree to a `multiPhylo` object whose own attributes apply to all trees,
#' for example trees read from a Nexus file, causes data to be lost.
#'
#' @param add Tree to add to the list, of class \code{\link[ape:read.tree]{phylo}}.
#' @param treeList A list of trees, of class \code{list},
#' \code{\link[ape:multiphylo]{multiPhylo}},
#' or, if a single tree, \code{\link[ape:read.tree]{phylo}}.
#'
#' @return A list of class \code{list} or \code{multiPhylo} (following the
#' original class of \code{treeList}), whose first element is the tree specified
#' as \code{add}.
#'
#' @template MRS
#'
#' @export
UnshiftTree <- function(add, treeList) {
  if (class(treeList) == 'multiPhylo') {
    structure(c(list(add), lapply(treeList, function (X) X)), class= 'multiPhylo')
  } else if (class(treeList) == 'phylo') {
    treeList <- structure(list(add, treeList), class='multiPhylo')
  } else { # including: if (class(trees) == 'list') {
    c(list(add), treeList)
  }
}
