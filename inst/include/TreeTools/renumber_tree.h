#ifndef TreeTools_renumber_tree_
#define TreeTools_renumber_tree_

#include <cassert>
#include <cstdio>
#include <cstdlib> /* for calloc */
#include <Rcpp.h>
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
      if (!parent_of[i]) root_node = i;
      if (!n_children[i]) ++n_tip;
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

    return (ret);
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

  // "Arkorder" is my term for a specific subset of postorder in which
  // edges are ordered such that all occurrences of each parent node
  // occur together.
  // Subtract one from $edge before passing.
  // [[Rcpp::export]]
  inline Rcpp::IntegerMatrix postorder_edges(
      const Rcpp::IntegerMatrix edge,
      const Rcpp::LogicalVector size_sort)
  {
    const int32
      n_edge = edge.nrow(),
      node_limit = n_edge + 1;

    int32
      root_node = 0,
      n_tip = 0;

    // 6 * checks we've enough memory for all children_of arrays too.
    // 0.9999 leaves room for memory overhead: seems in practice to avoid
    // attempting a doomed call to calloc.
    if (long(6 * node_limit * sizeof(int32)) > 0.9999L * INTPTR_MAX) {
      throw std::length_error("Tree too large for postorder_edges. "            // # nocov
                              "Try running 64-bit R?");                         // # nocov
    }

    int32 * parent_of = (int32*) std::calloc(node_limit, sizeof(int32)),
          * n_children = (int32*) std::calloc(node_limit, sizeof(int32)),
          * subtree_size = (int32*) std::calloc(node_limit, sizeof(int32));
    int32 ** children_of = new int32*[node_limit];

    for (int32 i = n_edge; i--; ) {
      parent_of[CHILD(i)] = PARENT(i);
      n_children[PARENT(i)] += 1;
    }

    for (int32 i = node_limit; i--; ) {
      if (!parent_of[i]) root_node = i;
      if (!n_children[i]) ++n_tip;
      children_of[i] = new int32[n_children[i]];
    }
    std::free(parent_of);

    int32 * found_children = (int32*) std::calloc(node_limit, sizeof(int32));
    for (int32 i = n_edge; i--; ) {
      children_of[PARENT(i)][found_children[PARENT(i)]] = CHILD(i);
      found_children[PARENT(i)] += 1;
    }
    std::free(found_children);

    const int32 n_node = n_edge - n_tip + 1;

    for (int32 tip = n_tip; tip--; ) {
      subtree_size[tip] = 1;
    }
    get_subtree_size(root_node, subtree_size, n_children, children_of, n_edge);

    for (int32 node = n_tip; node != node_limit; node++) {
      insertion_sort_by_smallest(children_of[node], n_children[node],
                                 subtree_size);
    }
    int32 * node_order = (int32*) malloc(n_node * sizeof(int32));
    for (int32 i = n_node; i--; ) {
      node_order[i] = i + n_tip;
    }
    if (size_sort[0]) {
      timsort_by_smallest(node_order, n_node, subtree_size);
    }
    std::free(subtree_size);

    Rcpp::IntegerMatrix ret(n_edge, 2);
    int32 this_edge = 0;
    for (int32 i = 0; i != n_node; ++i) {
      const int32 this_parent = node_order[i];
      for (int32 j = 0; j != n_children[this_parent]; ++j) {
        ret(this_edge, 0) = this_parent + 1;
        ret(this_edge, 1) = children_of[this_parent][j] + 1;
        ++this_edge;
      }
    }
    std::free(n_children);
    for (int32 i = node_limit; i--; ) {
      delete[] children_of[i];
    }
    delete[] (children_of);
    std::free(node_order);

    return (ret);
  }
}

#endif
