preorder_edges_and_nodes <- function(parent, child) {
  .Call(`_TreeTools_preorder_edges_and_nodes`, parent, child)
}

preorder_weighted <- function(parent, child, edgeLen) {
  .Call(`_TreeTools_preorder_weighted`, parent, child, edgeLen)
}

postorder_edges <- function(edge, sizeSort) {
  .Call(`_TreeTools_postorder_edges`, edge, sizeSort)
}

postorder_weighted <- function(edge, weight, sizeSort) {
  .Call(`_TreeTools_postorder_weighted`, edge, weight, sizeSort)
}

root_binary <- function(edge, outgroup) {
  .Call(`_TreeTools_root_binary`, edge, outgroup)
}

root_on_node <- function(phy, outgroup) {
  .Call(`_TreeTools_root_on_node`, phy, outgroup)
}
