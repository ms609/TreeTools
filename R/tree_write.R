#' Write a phylogenetic tree in Newick format
#'
#' Creates a character string describing a phylogenetic tree in Newick format,
#' using R's internal tip numbering.
#'
#'
#' @param x Object to convert to Newick format.
#' See Usage section for supported classes.
#'
#' @return `as.Newick` returns a character string representing `tree` in Newick
#' format.
#'
#' @seealso
#' - [`RenumberTips`]
#'
#' - [`ape::write.tree`]
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
