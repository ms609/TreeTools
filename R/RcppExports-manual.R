preorder_edges_and_nodes <- function(parent, child) {
  .Call(`_TreeTools_preorder_edges_and_nodes`, parent, child)
}

postorder_edges <- function(edge) {
  .Call(`_TreeTools_postorder_edges`, edge)
}

root_binary <- function(edge, outgroup) {
  .Call(`_TreeTools_root_binary`, edge, outgroup)
}

root_on_node <- function(phy, outgroup) {
  .Call(`_TreeTools_root_on_node`, phy, outgroup)
}
