#' Quick sample
#'
#' Faster than inbuilt sample because it avoids some checks. 
#' 
#' @param x A vector to sample.
#' @param len (optionally) integer specifying length of `x`.
#' 
#' @return A vector of length one, randomly sampled from `x`.
#' 
#' @examples 
#' SampleOne(9:10)
#' SampleOne(letters[1:4])
#' 
#' @template MRS
#' @keywords internal
#' @export
SampleOne <- function (x, len = length(x)) x[sample.int(len, 1L, FALSE, NULL, FALSE)]

Assert <- function (statement) if (!statement) stop(deparse(statement), " is FALSE")

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
  if (inherits(treeList, 'multiPhylo')) {
    structure(c(list(add), lapply(treeList, function (X) X)), class= 'multiPhylo')
  } else if (inherits(treeList, 'phylo')) {
    treeList <- structure(list(add, treeList), class='multiPhylo')
  } else { # including: if (class(trees) == 'list') {
    c(list(add), treeList)
  }
}
