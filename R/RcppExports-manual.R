keep_tip <- function(edge, keep) {
  .Call(`_BigTreeTools_keep_tip`, edge, keep)
}

preorder_edges_and_nodes <- function(parent, child) {
  .Call(`_BigTreeTools_preorder_edges_and_nodes`, parent, child)
}

preorder_weighted <- function(parent, child, edgeLen) {
  .Call(`_BigTreeTools_preorder_weighted`, parent, child, edgeLen)
}

postorder_order <- function(edge) {
  .Call(`_BigTreeTools_postorder_order`, edge)
}

root_binary <- function(edge, outgroup) {
  .Call(`_BigTreeTools_root_binary`, edge, outgroup)
}

root_on_node <- function(phy, outgroup) {
  .Call(`_BigTreeTools_root_on_node`, phy, outgroup)
}
