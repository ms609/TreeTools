#include <Rcpp/Lightest>
#include <stdexcept> /* for errors */
using namespace Rcpp;

#define R_TO_C 1

void mark_children(
    const int r_node_id,
    const IntegerVector& parent,
    const IntegerVector& child,
    LogicalVector& ret) {
  
  for (int i = parent.size(); i--; ) {
    if (parent[i] == r_node_id) {
      ret[i] = true;
      mark_children(child[i], parent, child, ret);
    }
  }
}

// [[Rcpp::export]]
LogicalMatrix descendant_edges(
    const IntegerVector parent,
    const IntegerVector child,
    const IntegerVector postorder) {
  if (parent.size() != child.size()) {
    Rcpp::stop("Parent and child must be the same length");
  }
  
  const int 
    n_edge = parent.size(),
    root = Rcpp::algorithm::min(parent.begin(), parent.end()),
    n_tip = root - 1,
    n_node = n_edge + 1
  ;
  LogicalMatrix ret(n_node - n_tip, n_edge);

  for (int i = 0; i != n_edge; i++) {
    const int
      edge_i = postorder[i] - R_TO_C,
      parent_i = parent[edge_i],
      child_i = child[edge_i],
      offset = R_TO_C + n_tip
    ;
    
    ret(parent_i - offset, edge_i) = true;
    if (child_i > n_tip) {
      for (int j = n_edge; j--; ) {
        if (ret(child_i - offset, j)) {
          ret(parent_i - offset, j) = true;
        }
      }
    }
  }
  return ret;
}
