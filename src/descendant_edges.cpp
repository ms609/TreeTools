#include <Rcpp/Lightest>
#include <cstdint> /* for uint64_t */
#include <stdexcept> /* for errors */
using namespace Rcpp;

#define R_TO_C 1

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
  
  const int root   = Rcpp::algorithm::min(parent.begin(), parent.end());
  const int n_tip  = root - 1;
  const int n_node = n_edge + 1;
  const int n_row  = n_node - n_tip;
  
  // Word-based storage (bitset) for intermediate computation
  const int n_words = (n_edge + 63) / 64;
  std::vector<std::vector<uint64_t>> bitmat(n_row, std::vector<uint64_t>(n_words, 0ULL));
  
  const int offset = R_TO_C + n_tip;
  
  // Traverse in postorder
  for (int i = 0; i < n_edge; i++) {
    const int edge_i   = postorder[i] - R_TO_C;
    const int parent_i = parent[edge_i] - offset;
    const int child_i  = child[edge_i]  - offset;
    
    // Mark edge itself
    bitmat[parent_i][edge_i / 64] |= (1ULL << (edge_i % 64));
    
    // Union child's descendants (if internal node)
    if (child_i >= 0) {
      for (int w = 0; w < n_words; w++) {
        bitmat[parent_i][w] |= bitmat[child_i][w];
      }
    }
  }
  
  // Convert back to LogicalMatrix for R (column-major fill)
  LogicalMatrix ret(n_row, n_edge);
  int* ret_ptr = INTEGER(ret);  // direct pointer to memory (0 = FALSE, 1 = TRUE)
  
  for (int j = 0; j < n_edge; j++) {
    const int word = j / 64;
    const int bit  = j % 64;
    const uint64_t mask = (1ULL << bit);
    
    for (int r = 0; r < n_row; r++) {
      ret_ptr[r + j * n_row] = (bitmat[r][word] & mask) ? 1 : 0;
    }
  }
  
  
  return ret;
}

// [[Rcpp::export]]
LogicalVector descendant_edges_single(
    const IntegerVector parent,
    const IntegerVector child,
    const IntegerVector postorder,
    const int edge_index,          // 1-based edge index from R
    const bool include_self = true // mimic DescendantEdges behaviour
) {
  const int n_edge = parent.size();
  if (child.size() != n_edge) {
    Rcpp::stop("`parent` and `child` must be the same length");
  }
  if (postorder.size() != n_edge) {
    Rcpp::stop("`postorder` must list each edge once");
  }
  if (edge_index < 1 || edge_index > n_edge) {
    Rcpp::stop("`edge_index` out of range");
  }
  
  LogicalVector ret(n_edge, false);
  std::vector<int> stack;
  stack.reserve(n_edge);
  
  // Start from the chosen edge
  const int start_child = child[edge_index - 1];
  if (include_self) {
    ret[edge_index - 1] = true;
  }
  stack.push_back(start_child);
  
  const int root = Rcpp::algorithm::min(parent.begin(), parent.end());
  const int n_tip = root - 1;
  
  // Depth-first traversal
  while (!stack.empty()) {
    int node = stack.back();
    stack.pop_back();
    
    if (node > n_tip) { // internal node
      for (int e = 0; e < n_edge; e++) {
        if (parent[e] == node) {
          ret[e] = true;
          stack.push_back(child[e]);
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
