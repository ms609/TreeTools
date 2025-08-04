#ifndef TreeTools_renumber_tree_
#define TreeTools_renumber_tree_

#include <cstdio>
#include <cstdlib> /* for calloc */
#include <stdexcept> /* for errors */
#include <Rcpp/Lightest>
#include "assert.h" /* for ASSERT */
#include "types.h"

#define MIN(a, b) ((a) < (b)) ? (a) : (b);
#define PARENT(i) edge[(i)]
#define CHILD(i) edge[(i) + n_edge]

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
                                         int32* sort_by) {
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

  // Modified from: https://www.geeksforgeeks.org/timsort/
  inline void merge(int32* arr, int32 l, int32 m, int32 r,
                    int32* sort_by) {

    const int32
      left_len = m - l + 1,
      right_len = r - m
    ;

    auto
      left = new int32[left_len],
      right = new int32[right_len]
    ;
    for (int32 i = left_len; i--; ) {
      left[i] = arr[l + i];
    }
    for (int32 i = right_len; i--; ) {
      right[i] = arr[m + 1 + i];
    }

    int32
      i = 0,
      j = 0,
      k = l
    ;

    while (i != left_len && j != right_len) {
      if (sort_by[left[i]] <= sort_by[right[j]]) {
        arr[k] = left[i];
        i++;
      } else {
        arr[k] = right[j];
        j++;
      }
      k++;
    }

    // Copy remaining elements of left, if any
    while (i != left_len) {
      arr[k] = left[i];
      ++k;
      ++i;
    }
    delete[] left;

    // Copy remaining element of right, if any
    while (j != right_len) {
      arr[k] = right[j];
      ++k;
      ++j;
    }
    delete[] right;
  }

  inline void timsort_by_smallest(int32* arr, int32 arr_len, int32* sort_by) {
    int32 run_length = arr_len;
    while (run_length > 64) ++run_length /= 2;

    // Sort individual subarrays of size run_length
    for (int32 i = 0; i < arr_len; i += run_length) {
      tim_insertion_sort_by_smallest(arr + i,
                                 i + run_length > arr_len ?
                                   arr_len % run_length : run_length,
                                 sort_by);
    }

    // Merge sorted subarrays
    for (int32 size = run_length; size < arr_len; size *= 2) {
      for (int32 left = 0; left < arr_len; left += 2 * size) {
        int32 mid = left + size - 1;
        int32 right = MIN(left + (2 * size) - 1, arr_len - 1);
        if (mid < right) merge(arr, left, mid, right, sort_by);
      }
    }
  }

  inline void add_child_edges(const int32 node, const int32 node_label,
                            int32 const* const* children_of,
                            const int32 *n_children,
                            Rcpp::IntegerMatrix& final_edges,
                            int32 *next_edge, int32 *next_label) {

    for (int32 child = 0; child != n_children[node]; ++child) {

      final_edges(*next_edge, 0) = node_label;
      const int32 this_child = children_of[node][child];

      if (n_children[this_child]) {

        const int32 child_label = *next_label;
        *next_label += 1;

        final_edges(*next_edge, 1) = child_label;
        *next_edge += 1;

        add_child_edges(this_child, child_label, children_of, n_children,
                        final_edges, next_edge, next_label);

      } else {

        final_edges(*next_edge, 1) = this_child;
        *next_edge += 1;

      }
    }
  }

  inline void add_child_edges(const int32 node, const int32 node_label,
                            int32 const* const* children_of,
                            const int32 *n_children,
                            const double *wt_above,
                            Rcpp::IntegerMatrix& final_edges,
                            Rcpp::NumericVector& final_weight,
                            int32 *next_edge, int32 *next_label) {

    for (int32 child = 0; child != n_children[node]; ++child) {

      final_edges(*next_edge, 0) = node_label;
      
      const int32 this_child = children_of[node][child];
      final_weight[*next_edge] = wt_above[this_child];
      
      if (n_children[this_child]) {

        const int32 child_label = *next_label;
        *next_label += 1;

        final_edges(*next_edge, 1) = child_label;
        *next_edge += 1;

        add_child_edges(this_child, child_label, children_of, n_children,
                        wt_above,
                        final_edges, final_weight, next_edge, next_label);

      } else {

        final_edges(*next_edge, 1) = this_child;
        *next_edge += 1;

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
    const int32
      n_edge = int32(parent.length()),
      node_limit = n_edge + 2
    ;

    if (child.length() != n_edge) {
      Rcpp::stop("Length of parent and child must match");
    }

    int32
      next_edge = 0,
      root_node = n_edge * 2, /* Initialize with too-big value */
      n_tip = 0
    ;

    int32 * parent_of = (int32*) std::calloc(node_limit, sizeof(int32)),
          * n_children = (int32*) std::calloc(node_limit, sizeof(int32)),
          * smallest_desc = (int32*) std::calloc(node_limit, sizeof(int32));
    int32 ** children_of = new int32*[node_limit];

    for (int32 i = n_edge; i--; ) {
      parent_of[child[i]] = parent[i];
      n_children[parent[i]] += 1;
    }

    for (int32 i = 1; i != node_limit; i++) {
      if (!parent_of[i]) {
        root_node = i;
      }
      if (!n_children[i]) {
        ++n_tip;
      } else {
        children_of[i] = new int32[n_children[i]];
      }
    }

    for (int32 tip = 1; tip != n_tip + 1; ++tip) {
      smallest_desc[tip] = tip;
      int32 parent = parent_of[tip];
      while (!smallest_desc[parent]) {
        smallest_desc[parent] = tip;
        parent = parent_of[parent];
      }
    }
    std::free(parent_of);

    int32 * found_children = (int32*) std::calloc(node_limit, sizeof(int32));
    for (int32 i = n_edge; i--; ) {
      children_of[parent[i]][found_children[parent[i]]] = child[i];
      found_children[parent[i]] += 1;
    }
    std::free(found_children);

    for (int32 node = n_tip + 1; node != node_limit; node++) {
      insertion_sort_by_smallest(children_of[node], n_children[node],
                                 smallest_desc);
    }
    std::free(smallest_desc);

    int32 next_label = n_tip + 2;
    Rcpp::IntegerMatrix ret(n_edge, 2);
    add_child_edges(root_node, n_tip + 1, children_of, n_children, ret,
                    &next_edge, &next_label);

    std::free(n_children);

    for (int32 i = n_tip + 1; i != node_limit; i++) {
      delete[] children_of[i];
    }
    delete[] children_of;

    return ret;
  }

// [[Rcpp::export]]
  inline Rcpp::List preorder_weighted(
      const Rcpp::IntegerVector parent,
      const Rcpp::IntegerVector child,
      const Rcpp::DoubleVector weight)
  {
    if (R_xlen_t(2LL + child.length() + 2LL + child.length()) >
          R_xlen_t(INT_FAST32_MAX)) {
      Rcpp::stop("Too many edges in tree: "                        // #nocov
                 "Contact 'TreeTools' maintainer for support.");   // #nocov
    }
    
    ASSERT(parent.length() < INT_FAST32_MAX - 2);
    const int32
      n_edge = int32(parent.length()),
      node_limit = n_edge + 2
    ;
    
    if (child.length() != n_edge) {
      Rcpp::stop("Length of parent and child must match");
    }
    if (weight.length() != n_edge) {
      Rcpp::stop("weights must match number of edges");
    }
    
    int32
      next_edge = 0,
      root_node = n_edge * 2, /* Initialize with too-big value */
      n_tip = 0
    ;
    
    int32 
      * parent_of = (int32*) std::calloc(node_limit, sizeof(int32)),
      * n_children = (int32*) std::calloc(node_limit, sizeof(int32)),
      * smallest_desc = (int32*) std::calloc(node_limit, sizeof(int32))
    ;
    double * wt_above = (double*) std::calloc(node_limit, sizeof(double));
    int32 ** children_of = new int32*[node_limit];
      
      for (int32 i = n_edge; i--; ) {
        wt_above[child[i]] = weight[i];
        parent_of[child[i]] = parent[i];
        n_children[parent[i]] += 1;
      }
      
      for (int32 i = 1; i != node_limit; i++) {
        if (!parent_of[i]) {
          root_node = i;
        }
        if (!n_children[i]) {
          ++n_tip;
        }
        children_of[i] = new int32[n_children[i]];
      }
      
      for (int32 tip = 1; tip != n_tip + 1; ++tip) {
        smallest_desc[tip] = tip;
        int32 parent = parent_of[tip];
        while (!smallest_desc[parent]) {
          smallest_desc[parent] = tip;
          parent = parent_of[parent];
        }
      }
      std::free(parent_of);
      
      int32 * found_children = (int32*) std::calloc(node_limit, sizeof(int32));
      for (int32 i = n_edge; i--; ) {
        children_of[parent[i]][found_children[parent[i]]] = child[i];
        found_children[parent[i]] += 1;
      }
      std::free(found_children);
      
      for (int32 node = n_tip + 1; node != node_limit; node++) {
        insertion_sort_by_smallest(children_of[node], n_children[node],
                                   smallest_desc);
      }
      std::free(smallest_desc);
      
      int32 next_label = n_tip + 2;
      Rcpp::IntegerMatrix ret(n_edge, 2);
      Rcpp::NumericVector ret_wt(n_edge);
      add_child_edges(root_node, n_tip + 1, children_of, n_children, wt_above,
                      ret, ret_wt, &next_edge, &next_label);
      
      std::free(wt_above);
      std::free(n_children);
      
      for (int32 i = 1; i != node_limit; i++) {
        delete[] children_of[i];
      }
      delete[] children_of;
      
      return Rcpp::List::create(ret, ret_wt);
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
    
    int32 stack_missing_children[STACK_THRESHOLD];
    bool stack_matched[STACK_THRESHOLD];
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
      missing_children = stack_missing_children;
      matched = stack_matched;
      ready_edges = stack_ready_edges;
      child_edges_start = stack_child_edges_start;
      child_edges = stack_child_edges;
      child_edge_counts = stack_child_edge_counts;
      
      std::memset(missing_children, 0, (node_limit + 1) * sizeof(int32));
      std::memset(matched, false, n_edge * sizeof(bool));
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
    for (int32 i = n_edge; i--; ) {
      ++missing_children[PARENT(i)];
      ++child_edge_counts[CHILD(i)];
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
        const int32 child = CHILD(i);
        child_edges[fill_pos[child]++] = i;
      }
    } else {
      int32* fill_pos = (int32*) std::malloc((node_limit + 1) * sizeof(int32));
      for (int32 node = 1; node <= node_limit; ++node) {
        fill_pos[node] = child_edges_start[node];
      }
      
      for (int32 i = 0; i < n_edge; ++i) {
        const int32 child = CHILD(i);
        child_edges[fill_pos[child]++] = i;
      }
      std::free(fill_pos);
    }
    
    // Find initial ready edges (leaves)
    int32 ready_count = 0;
    for (int32 i = n_edge; i--; ) {
      if (!missing_children[CHILD(i)]) {
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
      
      const int32 parent = PARENT(edge_idx);
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
      std::free(missing_children);
      std::free(matched);
      std::free(ready_edges);
      std::free(child_edges_start);
      std::free(child_edges);
      std::free(child_edge_counts);
    }
    
    return ret;
  }
}
#endif
