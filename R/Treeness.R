#' Relative length of internal branches
#' 
#' Treeness (also termed stemminess) is the proportion of a tree's length found
#' on internal branches \insertCite{Lanyon1988}{TreeTools}.
#' Insofar as external branches do not contain phylogenetic (grouping) signal,
#' trees with a high treeness can be interpreted as containing a higher
#' signal:noise ratio \insertCite{Phillips2003}{TreeTools}.
#' 
#' @template tree(s)Param
#' @returns `Treeness()` returns a numeric vector reporting the treeness of each
#' `tree`.
#' 
#' @examples
#' lowTree <- BalancedTree(6, lengths = c(1, 1, 4, 4, 4, 1, 1, 4, 4, 4))
#' plot(lowTree)
#' Treeness(lowTree)
#' highTree <- BalancedTree(6, lengths = c(6, 6, 1, 1, 1, 6, 6, 1, 1, 1))
#' plot(highTree)
#' Treeness(c(lowTree, highTree))
#' @template MRS
#' @family tree properties
#' @references \insertAllCited{}
#' @export
Treeness <- function(tree) {
  UseMethod("Treeness")
}

#' @export
Treeness.phylo <- function(tree) {
  if (is.null(tree[["edge.length"]])) {
    warning("No edge lengths specified")
    return(NULL)
  }
  edge <- tree[["edge"]]
  external <- edge[, 2] <= NTip(tree)
  weights <- tree[["edge.length"]]
  sum(weights[!external]) / sum(weights)
}

#' @export
Treeness.multiPhylo <- function(tree) {
  vapply(tree, Treeness.phylo, double(1))
}

#' @export
Treeness.list <- function(tree) {
  if (all(vapply(tree, inherits, TRUE, "phylo"))) {
    Treeness.multiPhylo(tree)
  } else {
    stop("`tree` must be a list of 'phylo' objects")
  }
}

#' @export
Treeness.NULL <- function(tree) NULL

#' @rdname Treeness
#' @export
Stemminess <- Treeness