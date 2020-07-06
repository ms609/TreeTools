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
SampleOne <- function (x, len = length(x)) {
  x[sample.int(len, 1L, FALSE, NULL, FALSE)]
}

#' Add tree to start of list
#'
#' `UnshiftTree()` adds a phylogenetic tree to the start of a list of trees.
#' This is useful where the class of a list of trees is unknown.
#'
#' Caution: adding a tree to a `multiPhylo` object whose own attributes apply
#' to all trees, for example trees read from a Nexus file, causes data to be
#' lost.
#'
#' @param add Tree to add to the list, of class [`phylo`][ape::read.tree].
#' @param treeList A list of trees, of class `list`,
#' [`multiPhylo`][ape::multiphylo], or, if a single tree,
#' [`phylo`][ape::read.tree].
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
UnshiftTree <- function (add, treeList) {
  if (inherits(treeList, 'multiPhylo')) {
    structure(c(list(add), lapply(treeList, I)), class = 'multiPhylo')
  } else if (inherits(treeList, 'phylo')) {
    structure(list(add, treeList), class = 'multiPhylo')
  } else { # including: if (is.list(trees)) {
    c(list(add), treeList)
  }
}

#' Apply a function that returns 64-bit integers over a list or vector
#'
#' Wrappers for members of the [`lapply()`] family intended for use when a
#' function `FUN` returns a vector of `integer64` objects.
#' `vapply()`, `sapply()` or `replicate()` drop the `integer64` class,
#' resulting in a vector of numerics that require conversion back to
#' 64-bit integers.  These functions restore the missing `class` attribute.
#'
#' @inheritParams base::lapply
#' @param FUN.LEN Integer specifying the length of the output of `FUN`.
#' @details For details of the underlying functions, see [`base::lapply()`].
#' @examples
#' sapply64(as.phylo(1:6, 6), as.TreeNumber)
#' vapply64(as.phylo(1:6, 6), as.TreeNumber, 1)
#' set.seed(0)
#' replicate64(6, as.TreeNumber(RandomTree(6)))
#' @template MRS
#' @seealso [`bit64::integer64()`][bit64-package]
#' @export
sapply64 <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  structure(sapply(X, FUN, ..., simplify, USE.NAMES), class = 'integer64')
}

#' @rdname sapply64
#' @export
vapply64 <- function (X, FUN, FUN.LEN = 1, ...) {
  structure(vapply(X, FUN, FUN.VALUE = numeric(FUN.LEN), ...),
            class = 'integer64')
}

#' @rdname sapply64
#' @export
replicate64 <- function (n, expr, simplify = "array") {
  sapply64(integer(n), eval.parent(substitute(function (...) expr)),
           simplify = simplify)
}
