#' Select element at random
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
#' @family utility functions
#' @export
SampleOne <- function(x, len = length(x)) {
  x[[sample.int(len, 1L, FALSE, NULL, FALSE)]]
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
#' @family utility functions
#'
#' @export
UnshiftTree <- function(add, treeList) {
  if (inherits(treeList, "multiPhylo")) {
    structure(c(list(add), lapply(treeList, I)), class = "multiPhylo")
  } else if (inherits(treeList, "phylo")) {
    structure(list(add, treeList), class = "multiPhylo")
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
#' @seealso \code{\link[bit64]{integer64}()}
#' @family utility functions
#' @export
sapply64 <- function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  structure(sapply(X, FUN, ..., simplify, USE.NAMES), class = "integer64")
}

#' @rdname sapply64
#' @export
vapply64 <- function(X, FUN, FUN.LEN = 1, ...) {
  structure(vapply(X, FUN, FUN.VALUE = numeric(FUN.LEN), ...),
            class = "integer64")
}

#' @rdname sapply64
#' @export
replicate64 <- function(n, expr, simplify = "array") {
  sapply64(integer(n), eval.parent(substitute(function(...) expr)),
           simplify = simplify)
}

#nocov start
#' Produce a legend for continuous gradient scales
#' 
#' Prints an annotated vertical bar coloured according to a continuous palette.
#' 
#' This function is now deprecated; it has been superseded by the more capable
#' [`PlotTools::SpectrumLegend()`] and will be removed in a future release.
# Deprecation notice added in TreeTools 1.9.2 (2023-04-25)
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
#' @param cex Character expansion factor relative to current `par("cex")`.
#' @param text.col Colour used for the legend text.
#' @param font,text.font Font used for the legend text; see [`text()`].
#' @param title Text to display 
#' @param title.col Colour for title; defaults to `text.col[1]`.
#' @param title.cex Expansion factor(s) for the title, defaults to `cex[1]`.
#' @param title.adj Horizontal adjustment for title: see the help for
#' `par("adj")`.
#' @param title.font Font used for the legend title.
#' @param pos,\dots Additional parameters to [`text()`].
#' 
#' @template MRS
#' @importFrom graphics segments strheight strwidth text
#' @keywords internal
#' @export
SpectrumLegend <- function(x0 = 0.05, y0 = 0.05,
                           x1 = x0, y1 = y0 + 0.2,
                           absolute = FALSE,
                           legend = character(0), palette,
                           lwd = 4, lty = 1, lend = "square", cex = 1,
                           text.col = par("col"),
                           font = NULL, text.font = font,
                           title = NULL, title.col = text.col[1], 
                           title.cex = cex[1], title.adj = 0.5, title.font = 2,
                           pos = 4,
                           ...) {
  
  .Deprecated("PlotTools::SpectrumLegend", package = "PlotTools")
  
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
       col = text.col,
       cex = cex,
       font = text.font,
       legend, pos = pos, ...)
  if (!is.null(title)) {
    text(mean(x0, x1) + (max(strwidth(legend)) / ifelse(pos == 2, -2, 2)),
         max(y0, y1) + prod(
           par("lheight"),
           strheight("")
         ),
         title,
         pos = 3,
         cex = title.cex, adj = title.adj, font = title.font, col = title.col,
         ...)
  }
}
#nocov end
