#ifndef TreeTools_renumber_tree_
#define TreeTools_renumber_tree_

#include <cstdio>
#include <cstdlib> /* for calloc */
#include <memory> /* for unique_ptr */
#include <stack> /* for stack */
#include <stdexcept> /* for errors */
#include <utility> /* for pair */
#include <vector>
#include <Rcpp/Lightest>
#include "assert.h" /* for ASSERT */
#include "types.h"

namespace TreeTools {
inline void swap(int32 *a, int32 *b) {
  const int32 temp = *a;
  *a = *b;
  *b = temp;
}

inline void insertion_sort_by_smallest(int32* arr, const int32 arr_len,
                                       const int32* sort_by) {
  ASSERT(arr_len > 0);
  switch (arr_len) {
  // case 0: return;
  case 1: return;
  case 2:
    if (sort_by[arr[0]] > sort_by[arr[1]]) {
      swap(&arr[0], &arr[1]);
    }
    return;
  }
  
  for (int32 i = 1; i != arr_len; ++i) {
    const int32 tmp = arr[i];
    const int32 key = sort_by[tmp];
    int32 j = i;
    while (j && sort_by[arr[j - 1]] > key) {
      arr[j] = arr[j - 1];
      --j;
    }
    arr[j] = tmp;
  }
}

// We'll use this a sentinel type to handle the unweighted case
struct NoWeights {};
// Used to conditionally create a type
struct DummyDoubleVector {};

struct TreeData {
  int32 n_edge;
  int32 node_limit;
  
  std::vector<int32> memory_block;
  
  int32* parent_of;
  int32* n_children;
  int32* smallest_desc;
  int32* children_start_idx;
  int32* children_data;
  
  TreeData(int32 num_edges)
    : n_edge(num_edges),
      node_limit(num_edges + 2),
      memory_block(
        node_limit + // parent_of
          node_limit + // n_children
          node_limit + // smallest_desc
          node_limit + // children_start_idx
          n_edge,      // children_data
          0)
  {
    auto it = memory_block.begin();
    parent_of = &*it; it += node_limit;
    n_children = &*it; it += node_limit;
    smallest_desc = &*it; it += node_limit;
    children_start_idx = &*it; it += node_limit;
    children_data = &*it;
  }
};

struct Frame {
  int32 node;
  int32 parent_label;
  int32 child_index;
  int32 child_count;
  const int32* node_children;
};

template <typename W, typename RetType>
RetType preorder_core(
    const Rcpp::IntegerVector& parent,
    const Rcpp::IntegerVector& child,
    const W& weights)
{
  const int32 n_edge = parent.length();
  if (R_xlen_t(2LL + child.length() + 2LL + child.length()) > R_xlen_t(INT_FAST32_MAX)) {
    Rcpp::stop("Too many edges in tree: Contact 'TreeTools' maintainer for support.");
  }
  
  if (child.length() != n_edge) {
    Rcpp::stop("Length of parent and child must match");
  }
  
  TreeData data(n_edge);
  int32 root_node = n_edge * 2;
  int32 n_tip = 0;
  
  std::conditional_t<std::is_same_v<W, NoWeights>, DummyDoubleVector, std::vector<double>> wt_above_storage;
  const std::vector<double>* wt_above_ptr = nullptr;
  
  if constexpr (!std::is_same_v<W, NoWeights>) {
    if (weights.length() != n_edge) {
      Rcpp::stop("weights must match number of edges");
    }
    wt_above_storage.resize(data.node_limit);
    wt_above_ptr = &wt_above_storage;
  }
  
  for (int32 i = 0; i < n_edge; ++i) {
    const int32 child_i = child[i];
    const int32 parent_i = parent[i];
    data.parent_of[child_i] = parent_i;
    ++data.n_children[parent_i];
    if constexpr (!std::is_same_v<W, NoWeights>) {
      wt_above_storage[child_i] = weights[i];
    }
  }
  
  int32 current_idx = 0;
  for (int32 i = 1; i < data.node_limit; i++) {
    if (!data.parent_of[i]) {
      root_node = i;
    }
    if (!data.n_children[i]) {
      ++n_tip;
    } else {
      data.children_start_idx[i] = current_idx;
      current_idx += data.n_children[i];
    }
  }
  
  for (int32 tip = 1; tip < n_tip + 1; ++tip) {
    data.smallest_desc[tip] = tip;
    int32 parent_node = data.parent_of[tip];
    while (parent_node && !data.smallest_desc[parent_node]) {
      data.smallest_desc[parent_node] = tip;
      parent_node = data.parent_of[parent_node];
    }
  }
  
  std::fill(data.n_children, data.n_children + data.node_limit, 0);
  for (int32 i = 0; i < n_edge; ++i) {
    int32 p = parent[i];
    int32 insert_pos = data.children_start_idx[p] + data.n_children[p];
    data.children_data[insert_pos] = child[i];
    ++data.n_children[p];
  }

  for (int32 node = n_tip + 1; node < data.node_limit; ++node) {
    int32* node_children = data.children_data + data.children_start_idx[node];
    insertion_sort_by_smallest(node_children, data.n_children[node],
                               data.smallest_desc);
  }
  
  int32 next_edge = 0;
  int32 next_label = n_tip + 2;
  
  Rcpp::IntegerMatrix ret_edges(n_edge, 2);
  
  std::conditional_t<std::is_same_v<W, NoWeights>, DummyDoubleVector, Rcpp::NumericVector> ret_weights;
  
  if constexpr (!std::is_same_v<W, NoWeights>) {
    ret_weights = Rcpp::NumericVector(n_edge);
  }
  
  std::stack<Frame> stack;
  int32 root_label = n_tip + 1;
  
  // Initialize with root node children
  {
    int32 child_count = data.n_children[root_node];
    if (child_count > 0) {
      stack.push(Frame{root_node, root_label, 0, child_count, data.children_data + data.children_start_idx[root_node]});
    }
  }
  
  while (!stack.empty()) {
    Frame& top = stack.top();
    
    if (top.child_index == top.child_count) {
      stack.pop();
      continue;
    }
    
    int32 child_node = top.node_children[top.child_index];
    
    ret_edges(next_edge, 0) = top.parent_label;
    if constexpr (!std::is_same_v<W, NoWeights>) {
      ret_weights[next_edge] = (*wt_above_ptr)[child_node];
    }
    
    if (data.n_children[child_node] == 0) {
      ret_edges(next_edge, 1) = child_node;
      ++next_edge;
      ++top.child_index;
    } else {
      int32 child_label = next_label++;
      ret_edges(next_edge, 1) = child_label;
      ++next_edge;
      ++top.child_index;
      
      int32 child_count = data.n_children[child_node];
      const int32* child_children = data.children_data + data.children_start_idx[child_node];
      stack.push(Frame{child_node, child_label, 0, child_count, child_children});
    }
  }
  
  if constexpr (std::is_same_v<RetType, Rcpp::IntegerMatrix>) {
    return ret_edges;
  } else {
    return std::make_pair(ret_edges, ret_weights);
  }
}
// === PUBLIC EXPORTED FUNCTIONS ===

// [[Rcpp::export]]
inline Rcpp::IntegerMatrix preorder_edges_and_nodes(
    const Rcpp::IntegerVector parent,
    const Rcpp::IntegerVector child)
{
  return preorder_core<NoWeights, Rcpp::IntegerMatrix>(parent, child, NoWeights{});
}

// [[Rcpp::export]]
inline Rcpp::List preorder_weighted(
    const Rcpp::IntegerVector& parent,
    const Rcpp::IntegerVector& child,
    const Rcpp::DoubleVector& weight)
{
  std::pair<Rcpp::IntegerMatrix, Rcpp::NumericVector> result =
    preorder_core<Rcpp::DoubleVector, std::pair<Rcpp::IntegerMatrix, Rcpp::NumericVector>>(parent, child, weight);
  
  Rcpp::List ret = Rcpp::List::create(
    Rcpp::Named("edge") = result.first,
    Rcpp::Named("edge.length") = result.second
  );
  
  return ret;
}








  inline int32 get_subtree_size(int32 node, int32 *subtree_size,
                                int32 *n_children, int32 **children_of,
                                int32 n_edge) {
    if (!subtree_size[node]) {
      for (int32 i = n_children[node]; i--; ) {
        subtree_size[node] += get_subtree_size(children_of[node][i],
                                subtree_size, n_children, children_of, n_edge);
      }
    }
    return subtree_size[node];
  }











  
  template <typename T, std::size_t StackSize>
  struct SmallBuffer {
    static_assert(std::is_trivial<T>::value,
                  "SmallBuffer requires a trivial type T");
    
    bool use_stack;
    std::array<T, StackSize> stack;
    T* heap;

    SmallBuffer(std::size_t needed)
      : use_stack(needed <= StackSize), heap(nullptr)
    {
      if (use_stack) {
        // zero only the number of elements we'll actually use
        std::memset(stack.data(), 0, needed * sizeof(T));
      } else {
        heap = static_cast<T*>(std::calloc(needed, sizeof(T)));
        if (!heap) throw std::bad_alloc{};
      }
    }
    
    ~SmallBuffer() {
      if (!use_stack) std::free(heap);
    }
    
    inline T* data() noexcept {
      return use_stack ? stack.data() : heap;
    }
  };
  
  // [[Rcpp::export]]
  inline Rcpp::IntegerVector postorder_order(const Rcpp::IntegerMatrix edge)
  {
    const int32 n_edge = edge.nrow();
    const int32 node_limit = n_edge + 1;
    
    if (long(6 * node_limit * sizeof(int32)) > 0.9999L * INTPTR_MAX) {
      Rcpp::stop("Tree too large for postorder_order. Try running 64-bit R?");
    }
    
    constexpr int32 STACK_THRESHOLD = 2048;
    SmallBuffer<int32, STACK_THRESHOLD> missing_children(node_limit + 1);
    SmallBuffer<char, STACK_THRESHOLD> matched(n_edge);
    
    // Count children
    int32* mc = missing_children.data();
    char*   m = matched.data();
    for (int32 i = 0; i < n_edge; ++i) {
      ++mc[edge[i]];
    }
    
    int32 found = 0;
    Rcpp::IntegerVector ret(n_edge);
    
    do {
      for (int32 i = n_edge; i--;) {
        if (!m[i] && !mc[edge[i + n_edge]]) {
          m[i] = true;
          --mc[edge[i]];
          ret[found++] = i + 1;
        }
      }
    } while (found != n_edge);
    
    return ret;
  }
}

#endif
