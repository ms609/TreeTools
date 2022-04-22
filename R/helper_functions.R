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
SampleOne <- function(x, len = length(x)) {
  x[sample.int(len, 1L, FALSE, NULL, FALSE)]
}

#' Add tree to start of list
#'
#' `UnshiftTree()` adds a phylogenetic tree to the start of a list of trees.
#' This is useful where the class of a list of trees is unknown, or where
#' names of trees should be retained.
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
#' @seealso
#' [`c()`] joins a tree or series of trees to a `multiPhylo` object, but loses
#' names and does not handle lists of trees.
#'
#' @examples
#' forest <- as.phylo(0:5, 6)
#' tree <- BalancedTree(6)
#'
#' UnshiftTree(tree, forest)
#' UnshiftTree(tree, tree)
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
sapply64 <- function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  structure(sapply(X, FUN, ..., simplify, USE.NAMES), class = 'integer64')
}

#' @rdname sapply64
#' @export
vapply64 <- function(X, FUN, FUN.LEN = 1, ...) {
  structure(vapply(X, FUN, FUN.VALUE = numeric(FUN.LEN), ...),
            class = 'integer64')
}

#' @rdname sapply64
#' @export
replicate64 <- function(n, expr, simplify = "array") {
  sapply64(integer(n), eval.parent(substitute(function(...) expr)),
           simplify = simplify)
}

#' Produce a legend for continuous gradient scales
#' 
#' Prints an annotated vertical bar coloured according to a continuous palette.
#' 
#' This convenience function is not yet very customizable; do file a GitHub
#' issue if you would value additional functionality.
#' 
#' @param x0,y0,x1,y1 Coordinates of the bottom-left and top-right end of the
#' bar.
#' @param absolute Logical specifying whether `x` and `y` values denote
#' coordinates (`TRUE`) or relative position, where (0, 0) denotes the
#' bottom-left of the plot area and (1, 1) the top right.
#' @param legend Character vector with which to label points on `palette`.
#' @param palette Colour palette to depict.
#' @param lwd,lty,lend Additional parameters to [`segments()`],
#' controlling line style.
#' @param pos,\dots Additional parameters to [`text()`].
#' 
#' @examples 
#' plot(0:1, 0:1, type = "n", frame.plot = FALSE,
#'      xlab = "x", ylab = "y")
#' SpectrumLegend(legend = c("Dark", "Middle", "Bright"),
#'                palette = hcl.colors(32L), lwd = 5)
#' SpectrumLegend(0.4, 0.95, 0.9, 0.95, abs = TRUE,
#'                legend = seq(1, 9, by = 2), palette = 1:9, pos = 1)
#' @template MRS
#' @importFrom graphics segments text
#' @export
SpectrumLegend <- function(x0 = 0.05, y0 = 0.05,
                           x1 = x0, y1 = y0 + 0.2,
                           absolute = FALSE,
                           legend = character(0), palette,
                           lwd = 4, lty = 1, lend = "square",
                           pos = 4, ...) {
  nCol <- length(palette)
  
  if (!absolute) {
    corners <- par("usr") # x0 x1 y0 y1
    xRange <- corners[2] - corners[1]
    yRange <- corners[4] - corners[3]
    
    # Order is important: lazy evaluation will set x1 = modified x0
    x1 <- corners[1] + (x1 * xRange)
    x0 <- corners[1] + (x0 * xRange)
    y1 <- corners[3] + (y1 * yRange)
    y0 <- corners[3] + (y0 * yRange)
  }
  
  segX <- x0 + ((x1 - x0) * 0:nCol / nCol)
  segY <- y0 + ((y1 - y0) * 0:nCol / nCol)
  
  nPlus1 <- nCol + 1L
  segments(segX[-nPlus1], segY[-nPlus1],
           segX[-1], segY[-1],
           col = palette,
           lwd = lwd, lty = lty, lend = lend)
  text(seq(x0, x1, length.out = length(legend)),
       seq(y0, y1, length.out = length(legend)),
       legend, pos = pos, ...)
}

