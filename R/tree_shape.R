#' Integer representing shape of a tree
#'
#' @template treeParam
#'
#' @return `TreeShape` returns an integer specifying the shape of a tree,
#' ignoring tip labels.
#'
#' @template MRS
#' @export
TreeShape <- function (tree) {
  tree <- root(tree, 1, resolve.root = TRUE)
  edge <- tree$edge
  nTip <- NTip(tree)
  edge <- PostorderEdges(edge[, 1], edge[, 2], nTip = nTip)

  edge_to_shape(edge[[1]], edge[[2]], nTip)

}
