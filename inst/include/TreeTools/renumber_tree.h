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

// Sentinel type to handle the unweighted case
struct NoWeights {};
// Used to conditionally create a type
struct DummyDoubleVector {};

struct TreeData {
  int32_t n_edge;
  int32_t node_limit;
  
  std::vector<int32_t> memory_block;
  
  int32_t* parent_of;
  int32_t* n_children;
  int32_t* smallest_desc;
  int32_t* children_start_idx;
  int32_t* children_data;
  
  TreeData(int32_t num_edges)
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
  int32_t node;
  int32_t parent_label;
  int32_t child_index;
  int32_t child_count;
  const int32_t* node_children;
};

struct PreorderState {
  TreeData& data;
  int32_t next_edge;
  int32_t next_label;
  int32_t n_tip;
  int32_t root_node;
  Rcpp::IntegerMatrix& ret_edges;
  
  PreorderState(TreeData& d, int32_t nt, int32_t rn, Rcpp::IntegerMatrix& edges)
    : data(d), next_edge(0), next_label(nt + 2), n_tip(nt), 
      root_node(rn), ret_edges(edges) {}
};

template<bool HasWeights>
inline void traverse_preorder(PreorderState& state, 
                              const double* wt_above = nullptr,
                              Rcpp::NumericVector* ret_weights = nullptr) {
  // Use a fixed-size stack for most trees to avoid heap allocation
  constexpr size_t STACK_SIZE = 128;
  std::array<Frame, STACK_SIZE> fast_stack;
  std::vector<Frame> heap_stack;
  
  size_t stack_pos = 0;
  bool use_heap = false;
  
  auto push_frame = [&](const Frame& f) {
    if (stack_pos < STACK_SIZE && !use_heap) {
      fast_stack[stack_pos++] = f;
    } else {
      if (!use_heap) {
        // Migrate to heap
        heap_stack.reserve(STACK_SIZE * 2);
        for (size_t i = 0; i < stack_pos; ++i) {
          heap_stack.push_back(fast_stack[i]);
        }
        use_heap = true;
      }
      heap_stack.push_back(f);
    }
  };
  
  auto pop_frame = [&]() -> Frame& {
    if (use_heap) {
      return heap_stack[heap_stack.size() - 1];
    } else {
      return fast_stack[stack_pos - 1];
    }
  };
  
  auto pop = [&]() {
    if (use_heap) {
      heap_stack.pop_back();
    } else {
      --stack_pos;
    }
  };
  
  auto empty = [&]() {
    return use_heap ? heap_stack.empty() : stack_pos == 0;
  };
  
  // Initialize with root
  int32_t child_count = state.data.n_children[state.root_node];
  if (child_count > 0) {
    push_frame({state.root_node, state.n_tip + 1, 0, child_count, 
               state.data.children_data + state.data.children_start_idx[state.root_node]});
  }
  
  while (!empty()) {
    Frame& top = pop_frame();
    
    if (top.child_index == top.child_count) {
      pop();
      continue;
    }
    
    int32_t child_node = top.node_children[top.child_index];
    
    state.ret_edges(state.next_edge, 0) = top.parent_label;
    
    if constexpr (HasWeights) {
      (*ret_weights)[state.next_edge] = wt_above[child_node];
    }
    
    if (state.data.n_children[child_node] == 0) {
      // Leaf node
      state.ret_edges(state.next_edge, 1) = child_node;
      ++state.next_edge;
      ++top.child_index;
    } else {
      // Internal node
      int32_t child_label = state.next_label++;
      state.ret_edges(state.next_edge, 1) = child_label;
      ++state.next_edge;
      ++top.child_index;
      
      push_frame({child_node, child_label, 0, state.data.n_children[child_node],
                 state.data.children_data + state.data.children_start_idx[child_node]});
    }
  }
}

// Separate functions to avoid template bloat in the setup code
inline Rcpp::IntegerMatrix preorder_unweighted_impl(
    const Rcpp::IntegerVector& parent,
    const Rcpp::IntegerVector& child) {
  
  const int32_t n_edge = parent.length();
  if (child.length() != n_edge) {
    Rcpp::stop("Length of parent and child must match");
  }
  
  TreeData data(n_edge);
  int32_t root_node = 0;
  int32_t n_tip = 0;
  
  for (int32_t i = n_edge; i--; ) {
    const int32_t child_i = child[i];
    const int32_t parent_i = parent[i];
    data.parent_of[child_i] = parent_i;
    ++data.n_children[parent_i];
  }
  
  int32_t current_idx = 0;
  for (int32_t i = 1; i < data.node_limit; i++) {
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
  
  for (int32_t tip = 1; tip < n_tip + 1; ++tip) {
    data.smallest_desc[tip] = tip;
    int32_t parent_node = data.parent_of[tip];
    while (parent_node && !data.smallest_desc[parent_node]) {
      data.smallest_desc[parent_node] = tip;
      parent_node = data.parent_of[parent_node];
    }
  }
  
  std::fill(data.n_children, data.n_children + data.node_limit, 0);
  for (int32_t i = 0; i < n_edge; ++i) {
    int32_t p = parent[i];
    int32_t insert_pos = data.children_start_idx[p] + data.n_children[p];
    data.children_data[insert_pos] = child[i];
    ++data.n_children[p];
  }
  
  for (int32_t node = n_tip + 1; node < data.node_limit; ++node) {
    int32_t* node_children = data.children_data + data.children_start_idx[node];
    std::sort(node_children, node_children + data.n_children[node],
              [&](int32_t a, int32_t b) {
                return data.smallest_desc[a] < data.smallest_desc[b];
              });
  }
  
  Rcpp::IntegerMatrix ret_edges(n_edge, 2);
  PreorderState state(data, n_tip, root_node, ret_edges);
  
  traverse_preorder<false>(state);
  
  return ret_edges;
}

inline Rcpp::List preorder_weighted_impl(
    const Rcpp::IntegerVector& parent,
    const Rcpp::IntegerVector& child,
    const Rcpp::DoubleVector& weights) {
  
  const int32_t n_edge = parent.length();
  if (child.length() != n_edge || weights.length() != n_edge) {
    Rcpp::stop("Length mismatch");
  }
  
  TreeData data(n_edge);
  std::vector<double> wt_above_storage(data.node_limit);
  int32_t root_node = 0;
  int32_t n_tip = 0;
  
  // Setup with weights
  for (int32_t i = n_edge; i--; ) {
    const int32_t child_i = child[i];
    const int32_t parent_i = parent[i];
    data.parent_of[child_i] = parent_i;
    ++data.n_children[parent_i];
    wt_above_storage[child_i] = weights[i];
  }
  
  // ... rest of setup same as unweighted ...
  int32_t current_idx = 0;
  for (int32_t i = 1; i < data.node_limit; i++) {
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
  
  for (int32_t tip = 1; tip < n_tip + 1; ++tip) {
    data.smallest_desc[tip] = tip;
    int32_t parent_node = data.parent_of[tip];
    while (parent_node && !data.smallest_desc[parent_node]) {
      data.smallest_desc[parent_node] = tip;
      parent_node = data.parent_of[parent_node];
    }
  }
  
  std::fill(data.n_children, data.n_children + data.node_limit, 0);
  for (int32_t i = 0; i < n_edge; ++i) {
    int32_t p = parent[i];
    int32_t insert_pos = data.children_start_idx[p] + data.n_children[p];
    data.children_data[insert_pos] = child[i];
    ++data.n_children[p];
  }
  
  for (int32_t node = n_tip + 1; node < data.node_limit; ++node) {
    int32_t* node_children = data.children_data + data.children_start_idx[node];
    std::sort(node_children, node_children + data.n_children[node],
              [&](int32_t a, int32_t b) {
                return data.smallest_desc[a] < data.smallest_desc[b];
              });
  }
  
  Rcpp::IntegerMatrix ret_edges(n_edge, 2);
  Rcpp::NumericVector ret_weights(n_edge);
  PreorderState state(data, n_tip, root_node, ret_edges);
  
  traverse_preorder<true>(state, wt_above_storage.data(), &ret_weights);
  
  return Rcpp::List::create(
    Rcpp::Named("edge") = ret_edges,
    Rcpp::Named("edge.length") = ret_weights
  );
}

// [[Rcpp::export]]
inline Rcpp::IntegerMatrix preorder_edges_and_nodes(
    const Rcpp::IntegerVector parent,
    const Rcpp::IntegerVector child) {
  return preorder_unweighted_impl(parent, child);
}

// [[Rcpp::export]]
inline Rcpp::List preorder_weighted(
    const Rcpp::IntegerVector& parent,
    const Rcpp::IntegerVector& child,
    const Rcpp::DoubleVector& weight) {
  return preorder_weighted_impl(parent, child, weight);
}

inline std::pair<Rcpp::IntegerMatrix, Rcpp::NumericVector> preorder_weighted_pair(
    const Rcpp::IntegerVector& parent,
    const Rcpp::IntegerVector& child,
    const Rcpp::DoubleVector& weights) {
  
  Rcpp::List result = preorder_weighted_impl(parent, child, weights);
  return std::make_pair(
    Rcpp::as<Rcpp::IntegerMatrix>(result["edge"]),
    Rcpp::as<Rcpp::NumericVector>(result["edge.length"])
  );
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
