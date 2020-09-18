preorder_edges_and_nodes <- function(parent, child) {
  .Call(`_TreeTools_preorder_edges_and_nodes`, parent, child)
}

postorder_edges <- function(edge) {
  .Call(`_TreeTools_postorder_edges`, edge)
}

root_on_node <- function(edge, outgroup) {
  .Call(`_TreeTools_root_on_node`, edge, outgroup)
}
