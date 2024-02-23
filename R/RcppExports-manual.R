keep_tip <- function(edge, keep) {
  .Call(`_TreeTools_keep_tip`, edge, keep)
}

preorder_edges_and_nodes <- function(parent, child) {
  .Call(`_TreeTools_preorder_edges_and_nodes`, parent, child)
}

preorder_weighted <- function(parent, child, edgeLen) {
  .Call(`_TreeTools_preorder_weighted`, parent, child, edgeLen)
}

postorder_order <- function(edge) {
  .Call(`_TreeTools_postorder_order`, edge)
}

root_binary <- function(edge, outgroup) {
  .Call(`_TreeTools_root_binary`, edge, outgroup)
}

#' Wrapper for internal C function `root_on_node()`
#' 
#' Direct entry point to `root_on_node()`; recommended for expert use only.
#' `RootTree()` checks that input is properly formatted and is recommended
#' for general use.
#' @param phy Minimally, a named list with entries `edge` and `Nnode`, in the
#' format of equivalent entries in a tree of class `phylo`. `edge.length` will
#' also be considered if supplied.
#' @param outgroup Numeric specifying index of leaf to use as outgroup.
#' @returns `root_on_node()` returns `phy` rooted on the speficied leaf.
#' 
#' @template MRS
#' @export
#' @keywords internal
root_on_node <- function(phy, outgroup) {
  .Call(`_TreeTools_root_on_node`, phy, outgroup)
}
