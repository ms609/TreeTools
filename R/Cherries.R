#' Count cherries in a tree
#' 
#' `Cherries()` counts the number of vertices in a binary tree whose children
#' are both leaves.
#' 
#' @param tree A binary tree, of class `phylo`; or a matrix corresponding to its
#' edge matrix.
#' @param nTip Number of leaves in tree.
#' @return `Cherries()` returns an integer specifying the number of nodes whose
#' children are both leaves.
#' @family tree properties
#' @template MRS
#' @export
Cherries <- function(tree, nTip) UseMethod("Cherries")

#' @rdname Cherries
#' @export
Cherries.phylo <- function(tree, nTip = NTip(tree)) {
  n_cherries_wrapper(tree[["edge"]][, 1], tree[["edge"]][, 2], nTip)
}

#' @rdname Cherries
#' @export
Cherries.numeric <- function(tree, nTip) {
  if (is.null(dim(tree)) || dim(tree)[[2]] != 2) {
    stop("`tree` must be the edge matrix of a tree of class phylo")
  }
  n_cherries_wrapper(tree[, 1], tree[, 2], nTip)
}
