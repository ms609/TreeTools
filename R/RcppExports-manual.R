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

#' C++ wrapper for root_on_node
#' 
#' Wrapper for C++ function `root_on_node()`.
#' No attempt is made to validate input; specifying invalid input may crash
#' R.
#' Behaviour is subject to change without notice.
#' 
#' @param phy Tree of class `phylo`. root node must == n_tip + 1.
#' @param outgroup Integer specifying new root node.
#' @template MRS
#' @export
root_on_node <- function(phy, outgroup) {
  .Call(`_TreeTools_root_on_node`, phy, outgroup)
}

#' C++ wrapper for keep_tip
#' 
#' Wrapper for C++ function `keep_tip()`.
#' No attempt is made to validate input; specifying invalid input may crash
#' R.
#' Behaviour is subject to change without notice.
#' 
#' @param edge Tree of class `phylo`, in (loose) preorder.
#' @param keep Logical vector specifying whether to keep each leaf.
#' @template MRS
#' @export
keep_tip <- function(edge, keep) {
  .Call(`_TreeTools_keep_tip`, edge, keep)
}
