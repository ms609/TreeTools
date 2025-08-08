#ifndef TreeTools_renumber_tree_
#define TreeTools_renumber_tree_

#include <cstdio>
#include <cstdlib> /* for calloc */
#include <stack> /* for stack */
#include <stdexcept> /* for errors */
#include <utility> /* for pair */
#include <Rcpp/Lightest>
#include "assert.h" /* for ASSERT */
#include "types.h"

namespace TreeTools {
  inline void swap(int32 *a, int32 *b) {
    const int32 temp = *a;
    *a = *b;
    *b = temp;
  }

  inline void quicksort_by_smallest(int32 *left, const int32 *right,
                                    const int32 *sort_by) {
    if (left >= right) return;

    const int32 pivot = sort_by[*right];
    int32 *centre = left;
    for (int32 *i = left; i <= right; i++) {
      if (sort_by[*i] <= pivot) {
        swap(centre, i);
        ++centre;
      }
    }
    quicksort_by_smallest(left, centre - 2, sort_by);
    quicksort_by_smallest(centre, right, sort_by);
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
      const int32
        tmp = arr[i],
        key = sort_by[tmp]
      ;
      int32 j = i;
      while (j && sort_by[arr[j - 1]] > key) {
        arr[j] = arr[j - 1];
        --j;
      }
      arr[j] = tmp;
    }
  }

  inline void insertion_sort_by_smallest(std::vector<int32>& arr,
                                         const int32 arr_len,
                                         const int32* sort_by) {
    ASSERT(arr_len > 0);
    switch (arr_len) {
    // case 0:
    case 1: return;
    case 2:
      if (sort_by[arr[0]] > sort_by[arr[1]]) {
        std::swap(arr[0], arr[1]);
      }
      return;
    default:
      for (size_t i = 1; i < size_t(arr_len); ++i) {
        const int32 tmp = arr[i];
        ASSERT(tmp >= 0 && tmp < ARR_LEN(sort_by));
        const int32 key = sort_by[tmp];
        size_t j = i;
        while (j > 0 && sort_by[arr[j - 1]] > key) {
          arr[j] = arr[j - 1];
          --j;
        }
        arr[j] = tmp;
      }
      
      // std::sort(arr.begin(), arr.end(), [&sort_by](int32 a, int32 b) {
      //   return sort_by[a] < sort_by[b];
      // });
    }
  }

  inline void tim_insertion_sort_by_smallest(int32* arr, const int32 arr_len,
                                             int32* sort_by) {
    ASSERT(arr_len > 0);
    for (int32 i = 1; i != arr_len; ++i) {
      const int32
        tmp = arr[i],
        key = sort_by[tmp]
      ;
      int32 j = i;
      while (j && sort_by[arr[j - 1]] > key) {
        arr[j] = arr[j - 1];
        --j;
      }
      arr[j] = tmp;
    }
  }

  struct Frame {
    int32 node;
    int32 parent_label;
    int32 child_index; // which child to process next
    int32 child_count;
    const int32* node_children;
  };

  inline void add_child_edges(
      const int32 root_node,
      const int32 root_label,
      const int32* children_data,
      const int32* children_start_idx,
      const int32* n_children,
      Rcpp::IntegerMatrix& ret,
      int32* next_edge,
      int32* next_label) {
    
    std::stack<Frame> stack;
    
    // Initialize with root node
    {
      int32 child_count = n_children[root_node];
      assert(child_count > 0);
      stack.push(Frame{root_node, root_label, 0, child_count, children_data + children_start_idx[root_node]});
    }
    
    while (!stack.empty()) {
      Frame& top = stack.top();
      
      if (top.child_index == top.child_count) {
        // All children processed for this node, pop stack
        stack.pop();
        continue;
      }
      
      int32 child = top.node_children[top.child_index];
      
      // Write edge info
      ret(*next_edge, 0) = top.parent_label;
      
      if (n_children[child] == 0) {
        ret(*next_edge, 1) = child;
        ++(*next_edge);
        ++top.child_index;
      } else {
        int32 child_label = *next_label;
        ret(*next_edge, 1) = child_label;
        ++(*next_label);
        ++(*next_edge);
        ++top.child_index;
        
        // Push child frame on stack to process its children next
        int32 child_count = n_children[child];
        const int32* child_children = children_data + children_start_idx[child];
        stack.push(Frame{child, child_label, 0, child_count, child_children});
      }
    }
  }

  inline void add_child_edges(
      const int32 root_node,
      const int32 root_label,
      const int32* children_data,
      const int32* children_start_idx,
      const int32* n_children,
      const std::vector<double>& wt_above,
      Rcpp::IntegerMatrix& ret,
      Rcpp::NumericVector& weight,
      int32* next_edge,
      int32* next_label) {
    
    std::stack<Frame> stack;
    
    // Initialize with root node
    {
      int32 child_count = n_children[root_node];
      assert(child_count > 0);
      stack.push(Frame{root_node, root_label, 0, child_count, children_data + children_start_idx[root_node]});
    }
    
    while (!stack.empty()) {
      Frame& top = stack.top();
      
      if (top.child_index == top.child_count) {
        // All children processed for this node, pop stack
        stack.pop();
        continue;
      }
      
      int32 child = top.node_children[top.child_index];
      
      // Write edge info
      ret(*next_edge, 0) = top.parent_label;
      weight[*next_edge] = wt_above[child];
      
      if (n_children[child] == 0) {
        ret(*next_edge, 1) = child;
        ++(*next_edge);
        ++top.child_index;
      } else {
        int32 child_label = *next_label;
        ret(*next_edge, 1) = child_label;
        ++(*next_label);
        ++(*next_edge);
        ++top.child_index;
        
        // Push child frame on stack to process its children next
        int32 child_count = n_children[child];
        const int32* child_children = children_data + children_start_idx[child];
        stack.push(Frame{child, child_label, 0, child_count, child_children});
      }
    }
  }
  
  // [[Rcpp::export]]
  inline Rcpp::IntegerMatrix preorder_edges_and_nodes(
      const Rcpp::IntegerVector parent,
      const Rcpp::IntegerVector child)
  {
    if (R_xlen_t(2LL + child.length() + 2LL + child.length()) > R_xlen_t(INT_FAST32_MAX)) {
      Rcpp::stop("Too many edges in tree: "                        // #nocov
                 "Contact 'TreeTools' maintainer for support.");   // #nocov
    }
    
    ASSERT(parent.length() < INT_FAST32_MAX - 2);
    const int32 n_edge = int32(parent.length());
    if (child.length() != n_edge) {
      Rcpp::stop("Length of parent and child must match");
    }
    const int32 max_node = n_edge + 1;
    assert(max_node == *std::max_element(parent.begin(), parent.end()));
    const int32 node_limit = max_node + 1;


    int32 next_edge = 0;
    int32 root_node = n_edge * 2; /* Initialize with too-big value */
    int32 n_tip = 0;
    
    // Single large allocation instead of many small ones
    const size_t total_ints_needed = 
      node_limit +                    // parent_of
      node_limit +                    // n_children  
      node_limit +                    // smallest_desc
      node_limit +                    // children_start_idx
      n_edge;                         // children_data
    
    std::vector<int32> memory_block(total_ints_needed, 0);
    
    auto it = memory_block.begin();
    int32* parent_of = &*it;
    it += node_limit;
    int32* n_children = &*it;
    it += node_limit;
    int32* smallest_desc = &*it;
    it += node_limit;
    int32* children_start_idx = &*it;
    it += node_limit;
    int32* children_data = &*it;
    
    for (int32 i = n_edge; i--; ) {
      const int32 parent_i = parent[i];
      parent_of[child[i]] = parent_i;
      ++n_children[parent_i];
    }
    
    int32 current_idx = 0;
    for (int32 i = 1; i < node_limit; i++) {
      if (!parent_of[i]) {
        root_node = i;
      }
      if (!n_children[i]) {
        ++n_tip;
      } else {
        children_start_idx[i] = current_idx;
        current_idx += n_children[i];
      }
    }

    for (int32 tip = 1; tip < n_tip + 1; ++tip) {
      smallest_desc[tip] = tip;
      int32 parent = parent_of[tip];
      while (!smallest_desc[parent]) {
        smallest_desc[parent] = tip;
        parent = parent_of[parent];
      }
    }
    
    // Reset n_children - use as insertion counter
    std::fill(n_children, n_children + node_limit, 0);
    for (int32 i = 0; i < n_edge; ++i) {
      int32 p = parent[i];
      int32 insert_pos = children_start_idx[p] + n_children[p];
      children_data[insert_pos] = child[i];
      ++n_children[p];
    }

    for (int32 node = n_tip + 1; node < node_limit; ++node) {
      int32* node_children = children_data + children_start_idx[node];
      insertion_sort_by_smallest(node_children, n_children[node], smallest_desc);
    }
    
    int32 next_label = n_tip + 2;
    Rcpp::IntegerMatrix ret(n_edge, 2);
    add_child_edges(root_node, n_tip + 1, children_data, children_start_idx,
                    n_children, ret, &next_edge, &next_label);

    return ret;
  }

  inline std::pair<Rcpp::IntegerMatrix, Rcpp::NumericVector> preorder_weighted_pair(
        const Rcpp::IntegerVector& parent,
        const Rcpp::IntegerVector& child,
        const Rcpp::DoubleVector& weight)
  {
    if (R_xlen_t(2LL + child.length() + 2LL + child.length()) >
          R_xlen_t(INT_FAST32_MAX)) {
      Rcpp::stop("Too many edges in tree: "                        // #nocov
                 "Contact 'TreeTools' maintainer for support.");   // #nocov
    }
    
    ASSERT(parent.length() < INT_FAST32_MAX - 2);
    const int32 n_edge = int32(parent.length());
    const int32 node_limit = n_edge + 2;
    
    if (child.length() != n_edge) {
      Rcpp::stop("Length of parent and child must match");
    }
    if (weight.length() != n_edge) {
      Rcpp::stop("weights must match number of edges");
    }
    
    int32 next_edge = 0;
    int32 root_node = n_edge * 2; /* Initialize with too-big value */
    int32 n_tip = 0;
    
    const size_t total_ints_needed = 
      node_limit +                    // parent_of
      node_limit +                    // n_children  
      node_limit +                    // smallest_desc
      node_limit +                    // children_start_idx
      n_edge;                         // children_data
    
    std::vector<int32> memory_block(total_ints_needed, 0);
    
    auto it = memory_block.begin();
    int32* parent_of = &*it;
    it += node_limit;
    int32* n_children = &*it;
    it += node_limit;
    int32* smallest_desc = &*it;
    it += node_limit;
    int32* children_start_idx = &*it;
    it += node_limit;
    int32* children_data = &*it;
    
    std::vector<double> wt_above(node_limit);
    
    for (int32 i = 0; i < n_edge; ++i) {
      const int32 child_i = child[i];
      const int32 parent_i = parent[i];
      wt_above[child_i] = weight[i];
      parent_of[child_i] = parent_i;
      ++n_children[parent_i];
    }
    
    int32 current_idx = 0;
    for (int32 i = 1; i < node_limit; i++) {
      if (!parent_of[i]) {
        root_node = i;
      }
      if (!n_children[i]) {
        ++n_tip;
      } else {
        children_start_idx[i] = current_idx;
        current_idx += n_children[i];
      }
    }
    
    for (int32 tip = 1; tip < n_tip + 1; ++tip) {
      smallest_desc[tip] = tip;
      int32 parent = parent_of[tip];
      while (!smallest_desc[parent]) {
        smallest_desc[parent] = tip;
        parent = parent_of[parent];
      }
    }
    
    // Reset n_children - use as insertion counter
    std::fill(n_children, n_children + node_limit, 0);
    for (int32 i = 0; i < n_edge; ++i) {
      int32 p = parent[i];
      int32 insert_pos = children_start_idx[p] + n_children[p];
      children_data[insert_pos] = child[i];
      ++n_children[p];
    }
    
    for (int32 node = n_tip + 1; node < node_limit; ++node) {
      int32* node_children = children_data + children_start_idx[node];
      insertion_sort_by_smallest(node_children, n_children[node], smallest_desc);
    }
    
    int32 next_label = n_tip + 2;
    Rcpp::IntegerMatrix ret(n_edge, 2);
    Rcpp::NumericVector ret_wt(n_edge);
    add_child_edges(root_node, n_tip + 1, children_data, children_start_idx,
                    n_children, wt_above, ret, ret_wt, &next_edge, &next_label);
    
    return std::make_pair(ret, ret_wt);
  }

  inline Rcpp::List preorder_weighted(
      const Rcpp::IntegerVector& parent,
      const Rcpp::IntegerVector& child,
      const Rcpp::DoubleVector& weight)
  {
    // Call the core function to get the pair
    std::pair<Rcpp::IntegerMatrix, Rcpp::NumericVector> result = 
      preorder_weighted_pair(parent, child, weight);
    
    // Manually create an Rcpp::List and populate it
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

  // [[Rcpp::export]]
  inline Rcpp::IntegerVector postorder_order(const Rcpp::IntegerMatrix edge)
  {
    const int32 n_edge = edge.nrow();
    const int32 node_limit = n_edge + 1;
    
    if (long(6 * node_limit * sizeof(int32)) > 0.9999L * INTPTR_MAX) {
      Rcpp::stop("Tree too large for postorder_order. "            // # nocov
                   "Try running 64-bit R?");                       // # nocov
    }
    
    constexpr int32 STACK_THRESHOLD = 2048;
    const bool use_stack = node_limit < STACK_THRESHOLD;
    std::array<int32, STACK_THRESHOLD> stack_missing_children;
    std::array<bool, STACK_THRESHOLD>  stack_matched;
    
    int32 stack_ready_edges[STACK_THRESHOLD];
    int32 stack_child_edges_start[STACK_THRESHOLD];
    int32 stack_child_edges[STACK_THRESHOLD];
    int32 stack_child_edge_counts[STACK_THRESHOLD];
    
    int32 * missing_children;
    bool* matched;
    int32* ready_edges;
    int32* child_edges_start;  // Start index for each node's child edges
    int32* child_edges;        // Flat array of child edges
    int32* child_edge_counts;  // Count of edges ending at each node
    
    if (use_stack) {
      missing_children = stack_missing_children.data();
      matched = stack_matched.data();
      ready_edges = stack_ready_edges;
      child_edges_start = stack_child_edges_start;
      child_edges = stack_child_edges;
      child_edge_counts = stack_child_edge_counts;
      
      std::fill(missing_children, missing_children + node_limit + 1, 0);
      std::fill(matched, matched + n_edge, false);
      std::memset(child_edges_start, 0, (node_limit + 1) * sizeof(int32));
      std::memset(child_edge_counts, 0, (node_limit + 1) * sizeof(int32));
    } else {
      missing_children = (int32*) std::calloc(node_limit + 1, sizeof(int32));
      matched = (bool*) std::calloc(n_edge, sizeof(bool));
      ready_edges = (int32*) std::malloc(n_edge * sizeof(int32));
      child_edges_start = (int32*) std::calloc(node_limit + 1, sizeof(int32));
      child_edges = (int32*) std::malloc(n_edge * sizeof(int32));
      child_edge_counts = (int32*) std::calloc(node_limit + 1, sizeof(int32));
    }
    
    // Count children for each node AND count edges ending at each node
    for (int32 i = 0; i < n_edge; ++i) {
      ++missing_children[edge[i]];
      ++child_edge_counts[edge[i + n_edge]];
    }
    
    // Build start indices for child_edges array
    int32 total = 0;
    for (int32 node = 1; node <= node_limit; ++node) {
      child_edges_start[node] = total;
      total += child_edge_counts[node];
    }
    
    // Populate child_edges array
    if (use_stack) {
      int32 fill_pos[STACK_THRESHOLD];
      for (int32 node = 1; node <= node_limit; ++node) {
        fill_pos[node] = child_edges_start[node];
      }
      
      for (int32 i = 0; i < n_edge; ++i) {
        const int32 child = edge[i];
        child_edges[fill_pos[child]++] = i;
      }
    } else {
      int32* fill_pos = (int32*) std::malloc((node_limit + 1) * sizeof(int32));
      for (int32 node = 1; node <= node_limit; ++node) {
        fill_pos[node] = child_edges_start[node];
      }
      
      for (int32 i = 0; i < n_edge; ++i) {
        const int32 child = edge[i];
        child_edges[fill_pos[child]++] = i;
      }
      std::free(fill_pos);
    }
    
    // Find initial ready edges (leaves)
    int32 ready_count = 0;
    for (int32 i = n_edge; i--; ) {
      if (!missing_children[edge[i]]) {
        ready_edges[ready_count++] = i;
      }
    }
    
    // Process ready edges
    int32 found = 0;
    int32 ready_pos = 0;
    Rcpp::IntegerVector ret(n_edge);

    while (ready_pos < ready_count) {
      const int32 edge_idx = ready_edges[ready_pos++];
      matched[edge_idx] = true;
      ret[found++] = edge_idx + 1;
      
      const int32 parent = edge[edge_idx + n_edge];
      if (--missing_children[parent] == 0) {
        // Parent node is now ready - add all edges ending at parent
        const int32 start = child_edges_start[parent];
        const int32 end = (parent == node_limit) ? n_edge : child_edges_start[parent + 1];
        
        for (int32 j = start; j < end; ++j) {
          const int32 candidate_edge = child_edges[j];
          if (!matched[candidate_edge]) {
            ready_edges[ready_count++] = candidate_edge;
          }
        }
      }
    }

    if (!use_stack) {
      std::free(ready_edges);
      std::free(child_edges_start);
      std::free(child_edges);
      std::free(child_edge_counts);
    }
    return ret;
  }
}
#endif
