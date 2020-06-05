#' Write a phylogenetic tree in Newick format
#'
#' `as.Newick()` creates a character string representation of a phylogenetic
#' tree, in the Newick format, using R's internal tip numbering.
#' Use [`RenumberTips()`] to ensure that the internal numbering follows the
#' order you expect.
#'
#'
#' @param x Object to convert to Newick format.
#' See Usage section for supported classes.
#'
#' @return `as.Newick()` returns a character string representing `tree` in Newick
#' format.
#'
#' @examples
#' trees <- list(BalancedTree(1:8), PectinateTree(8:1))
#' trees <- lapply(trees, RenumberTips, 1:8)
#' as.Newick(trees)
#'
#' @seealso
#' - Retain leaf labels: [`NewickTree()`]
#'
#' - Change R's internal numbering of leaves: [`RenumberTips()`]
#'
#' - Write tree to text or file: [`ape::write.tree()`]
#'
#' @template MRS
#' @export
as.Newick <- function (x) UseMethod('as.Newick')

#' @rdname as.Newick
#' @export
as.Newick.phylo <- function (x) {
  as_newick(Preorder(x)$edge - 1L)
}

#' @rdname as.Newick
#' @export
as.Newick.list <- function (x) {
  vapply(x, as.Newick, character(1L))
}

#' @rdname as.Newick
#' @export
as.Newick.multiPhylo <- as.Newick.list
