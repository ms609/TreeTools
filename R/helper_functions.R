#' Quickly sample
#'
#' `SampleOne()` is a fast alternative to  [`sample()`] that avoids some checks.
#'
#' @param x A vector to sample.
#' @param len (Optional) Integer specifying length of `x`.
#'
#' @return `SampleOne()` returns a length one vector, randomly sampled from `x`.
#'
#' @examples
#' SampleOne(9:10)
#' SampleOne(letters[1:4])
#'
#' @template MRS
#' @keywords internal
#' @export
SampleOne <- function (x, len = length(x)) x[sample.int(len, 1L, FALSE, NULL, FALSE)]

#' Add tree to start of list
#'
#' `UnshiftTree()` adds a phylogenetic tree to the start of a list of trees.
#' This is useful where the class of a list of trees is unknown.
#'
#' Caution: adding a tree to a `multiPhylo` object whose own attributes apply
#' to all trees, for example trees read from a Nexus file, causes data to be
#' lost.
#'
#' @param add Tree to add to the list, of class [`phylo`][ape:read.tree].
#' @param treeList A list of trees, of class `list`,
#' [`multiPhylo`][ape:multiphylo], or, if a single tree,
#' [`phylo`][ape:read.tree].
#'
#' @return `UnshiftTree()` returns a list of class `list` or `multiPhylo`
#' (following the original class of `treeList`), whose first element is the
#' tree specified as `add.
#'
#' @examples
#' forest <- as.phylo(0:5, 6)
#' tree <- BalancedTree(6)
#'
#' UnshiftTree(tree, forest)
#' UnshiftTree(tree, tree)
#'
#' @template MRS
#'
#' @export
UnshiftTree <- function(add, treeList) {
  if (inherits(treeList, 'multiPhylo')) {
    structure(c(list(add), lapply(treeList, I)), class = 'multiPhylo')
  } else if (inherits(treeList, 'phylo')) {
    structure(list(add, treeList), class = 'multiPhylo')
  } else { # including: if (is.list(trees)) {
    c(list(add), treeList)
  }
}
