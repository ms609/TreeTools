#include <Rcpp/Lightest>
#include <stdexcept> /* for errors */
using namespace Rcpp;

constexpr int R_TO_C = 1;

// [[Rcpp::export]]
LogicalMatrix descendant_edges(
    const IntegerVector parent,
    const IntegerVector child,
    const IntegerVector postorder) {
  
  const int n_edge = parent.size();
    
  if (child.size() != n_edge) {
    Rcpp::stop("`parent` and `child` must be the same length");
  }
  if (postorder.size() != n_edge) {
    Rcpp::stop("`postorder` must list each edge once");
  }
  
  const int 
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

// [[Rcpp::export]]
LogicalMatrix descendant_tips(
    const IntegerVector parent,
    const IntegerVector child,
    const IntegerVector postorder) {
  
  const int n_edge = parent.size();
    
  if (child.size() != n_edge) {
    Rcpp::stop("`parent` and `child` must be the same length");
  }
  if (postorder.size() != n_edge) {
    Rcpp::stop("`postorder` must list each edge once");
  }
  
  const int
    root = Rcpp::algorithm::min(parent.begin(), parent.end()),
    n_tip = root - 1,
    n_node = n_edge + 1
  ;
  LogicalMatrix ret(n_node, n_tip);

  for (int i = 0; i != n_edge; i++) {
    const int
      edge_i = postorder[i] - R_TO_C,
      parent_i = parent[edge_i],
      child_i = child[edge_i],
      offset = R_TO_C
    ;
    
    if (child_i > n_tip) {
      for (int j = n_tip; j--; ) {
        if (ret(child_i - offset, j)) {
          ret(parent_i - offset, j) = true;
        }
      }
    } else {
      ret(child_i - offset, child_i - R_TO_C) = true;
      ret(parent_i - offset, child_i - R_TO_C) = true;
    }
  }
  return ret;
}
