#ifndef TreeTools_renumber_tree_
#define TreeTools_renumber_tree_

#include <cassert>
#include <cstdio>
#include <cstdlib> /* for calloc */
#include <Rcpp/Lightest>
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
    assert(arr_len > 0);
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
    assert(arr_len > 0);
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
    if (2.0 * (2 + child.length()) > double(INT_FAST32_MAX)) {
      throw std::length_error("Too many edges in tree: "                        // #nocov
                              "Contact 'TreeTools' maintainer for support.");   // #nocov
    }

    const int32
      n_edge = parent.length(),
      node_limit = n_edge + 2
    ;

    if (child.length() != n_edge) {
      throw std::invalid_argument("Length of parent and child must match");
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
    add_child_edges(root_node, n_tip + 1, children_of, n_children, ret,
                    &next_edge, &next_label);

    std::free(n_children);

    for (int32 i = 1; i != node_limit; i++) {
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
    if (2.0 * (2 + child.length()) > double(INT_FAST32_MAX)) {
      throw std::length_error("Too many edges in tree: "                        // #nocov
                              "Contact 'TreeTools' maintainer for support.");   // #nocov
    }
    
    const int32
      n_edge = parent.length(),
      node_limit = n_edge + 2
    ;
    
    if (child.length() != n_edge) {
      throw std::invalid_argument("Length of parent and child must match");
    }
    if (weight.length() != n_edge) {
      throw std::invalid_argument("weights must match number of edges");
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
    const int32
      n_edge = edge.nrow(),
      node_limit = n_edge + 1;
    
    
    if (long(6 * node_limit * sizeof(int32)) > 0.9999L * INTPTR_MAX) {
      throw std::length_error("Tree too large for postorder_order. "            // # nocov
                              "Try running 64-bit R?");                         // # nocov
    }
    
    int32 * missing_children = (int32*) std::calloc(node_limit + 1, sizeof(int32));
    for (int32 i = n_edge; i--; ) {
      ++missing_children[PARENT(i)];
    }
    
    int32 found = 0;
    bool * matched = (bool*) std::calloc(node_limit, sizeof(bool));
    Rcpp::IntegerVector ret(n_edge);
    do {
      for (int32 i = n_edge; i--; ) {
        if (!matched[i]) {
          if (!missing_children[CHILD(i)]) {
            matched[i] = true;
            --missing_children[PARENT(i)];
            ret[found++] = i + 1;
          }
        }
      }
    } while (found != n_edge);
    std::free(missing_children);
    std::free(matched);
    
    return ret;
  }
}

#endif
