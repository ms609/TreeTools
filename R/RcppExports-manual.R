# Manual exports are required for functions that are defined in the 
# TreeTools namespace in header files.
keep_tip <- function(edge, keep) {
  .Call(`_TreeTools_keep_tip`, edge, keep)
}

#' Validate preorder edge ordering
#'
#' Checks whether edges are in valid preorder: each edge's parent has been
#' introduced (appeared as a child) by an earlier edge, or is the root.
#'
#' @param parent,child Integer vectors of edge endpoints (1-indexed).
#' @param nTip Integer: number of tips in the tree.
#' @returns Logical: `TRUE` if edges are in valid preorder, `FALSE` otherwise.
#' @examples
#' tree <- BalancedTree(8)
#' is_valid_preorder(tree[["edge"]][, 1], tree[["edge"]][, 2], NTip(tree))
#' @template MRS
#' @export
is_valid_preorder <- function(parent, child, nTip) {
  .Call(`_TreeTools_is_valid_preorder`, parent, child, nTip)
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
#' @param outgroup Integer specifying index of leaf or node to set as the
#' outgroup.
#' @returns `root_on_node()` returns `phy` rooted on the specified node.
#' 
#' @template MRS
#' @export
#' @keywords internal
root_on_node <- function(phy, outgroup) {
  .Call(`_TreeTools_root_on_node`, phy, outgroup)
}
