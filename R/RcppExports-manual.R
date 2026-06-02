# Manual exports are required for functions that are defined in the 
# TreeTools namespace in header files.
keep_tip <- function(edge, keep) {
  .Call(`_TreeTools_keep_tip`, edge, keep)
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
#'
#' `phy` may be supplied in any valid edge order; it is preordered internally,
#' so there is no need to call [`Preorder()`] first.  The returned tree is
#' always in preorder, with node numbers reassigned accordingly.
#' @param phy Minimally, a named list with entries `edge` and `Nnode`, in the
#' format of equivalent entries in a tree of class `phylo`. `edge.length` will
#' also be considered if supplied. The root node must be numbered `n_tip + 1`,
#' per the `phylo` convention.
#' @param outgroup Integer specifying index of leaf or node to set as the
#' outgroup.
#' @returns `root_on_node()` returns `phy` rooted on the specified node, in
#' preorder.
#' 
#' @template MRS
#' @export
#' @keywords internal
root_on_node <- function(phy, outgroup) {
  .Call(`_TreeTools_root_on_node`, phy, outgroup)
}
