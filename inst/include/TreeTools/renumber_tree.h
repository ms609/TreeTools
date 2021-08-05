#ifndef TreeTools_renumber_tree_
#define TreeTools_renumber_tree_

#include <cstdio>
#include <cstdlib> /* for calloc */
#include <Rcpp.h>
#include "types.h"
using namespace Rcpp;

namespace TreeTools {
  inline void swap(int32 *a, int32 *b) {
    const int32 temp = *a;
    *a = *b;
    *b = temp;
  }

  /* Requires unsigned integers. */
  /* If we chose signed, we'd have to impose a limit on n_children, which
   * would exclude star trees */
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

  inline void add_child_edges(const int32 node, const int32 node_label,
                            int32 const* const* children_of,
                            const int32 *n_children,
                            IntegerMatrix final_edges,
                            int32 *next_edge, int32 *next_label,
                            const int32 *n_tip, const int32 *n_edge) {

    for (int32 child = 0; child != n_children[node]; child++) {

      final_edges(*next_edge, 0) = node_label;
      int32 this_child = children_of[node][child];

      if (this_child <= *n_tip) {

        final_edges(*next_edge, 1) = this_child;
        ++(*next_edge);

      } else {

        const int32 child_label = (*next_label)++;
        final_edges(*next_edge, 1) = child_label;
        ++(*next_edge);

        add_child_edges(this_child, child_label, children_of, n_children,
                        final_edges, next_edge, next_label, n_tip, n_edge);

      }
    }
  }

  // [[Rcpp::export]]
  inline IntegerMatrix preorder_edges_and_nodes(const IntegerVector parent,
                                                const IntegerVector child)
  {
    if (2.0 * (2 + child.length()) > double(INT_FAST32_MAX)) {
      throw std::length_error("Too many edges in tree: "
                              "Contact 'TreeTools' maintainer for support.");
    }

    const int32 n_edge = parent.length(),
               node_limit = n_edge + 2;

    if (child.length() != n_edge) {
      throw std::invalid_argument("Length of parent and child must match");
    }

    int32 next_edge = 0,
         root_node = n_edge * 2, /* Initialize with too-big value */
         n_tip = 0;

    int32 * parent_of = (int32*) std::calloc(node_limit, sizeof(int32)),
          * n_children = (int32*) std::calloc(node_limit, sizeof(int32)),
          * smallest_desc = (int32*) std::calloc(node_limit, sizeof(int32));
    int32 ** children_of = new int32*[node_limit];

    for (int32 i = n_edge; i--; ) {
      parent_of[child[i]] = parent[i];
      ++(n_children[parent[i]]);
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
      children_of[parent[i]][(found_children[parent[i]])++] = child[i];
    }
    std::free(found_children);

    for (int32 node = n_tip + 1; node != node_limit; node++) {
      quicksort_by_smallest(children_of[node],
                            children_of[node] + n_children[node] - 1,
                            smallest_desc);
    }
    std::free(smallest_desc);

    int32 next_label = n_tip + 2;
    IntegerMatrix ret(n_edge, 2);
    add_child_edges(root_node, n_tip + 1,
                    children_of, n_children, ret,
                    &next_edge, &next_label, &n_tip, &n_edge);

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
  inline IntegerMatrix postorder_edges(const IntegerMatrix edge)
  {
    if (1L + edge.nrow() > long(0x7FFF)) {
      throw std::length_error("Too many edges in tree for postorder_edges: "
                              "Contact maintainer for advice");
      // In theory we could use INTX_MAX, which is larger than 16 bits on linux,
      // or we could change int32 to 32 bit.  The former has caused a seg fault
      // with invalid permissions on linux builds, possibly related to callocing
      // children_of?
      // Rather than attempt to debug now, I've chosen to place a hard limit on
      // edge.nrow for the time being.  --MS, 2020-05-26
      // Possibly fixed in branch optim-postorder, to check: 2021-08-05
    }

    const int32
      n_edge = edge.nrow(),
      node_limit = n_edge + 1;

    int32
      root_node = 0,
      n_tip = 0;

    // 6 * checks we've enough memory for all children_of arrays too.
    // 0.9999 leaves room for memory overhead: seems in practice to avoid
    // attempting a doomed call to calloc.
    if (long(6 * node_limit * sizeof(intx)) > 0.9999L * INTPTR_MAX) {
      throw std::length_error("Tree too large for postorder_edges. "            // # nocov
                              "Try running 64-bit R?");                         // # nocov
    }

    int32 * parent_of = (int32*) std::calloc(node_limit, sizeof(int32)),
          * n_children = (int32*) std::calloc(node_limit, sizeof(int32)),
          * subtree_size = (int32*) std::calloc(node_limit, sizeof(int32));
    int32 ** children_of = new int32*[node_limit];

    for (int32 i = n_edge; i--; ) {
      parent_of[edge(i, 1)] = edge(i, 0);
      ++(n_children[edge(i, 0)]);
    }

    for (int32 i = node_limit; i--; ) {
      if (parent_of[i] == 0) root_node = i;
      if (n_children[i] == 0) ++n_tip;
      children_of[i] = new int32[n_children[i]];
    }
    std::free(parent_of);

    int32 * found_children = (int32*) std::calloc(node_limit, sizeof(int32));
    for (int32 i = n_edge; i--; ) {
      children_of[edge(i, 0)][(found_children[edge(i, 0)])++] = edge(i, 1);
    }
    std::free(found_children);

    const int32 n_node = n_edge - n_tip + 1;

    for (int32 tip = 0; tip != n_tip; tip++) {
      subtree_size[tip] = 1;
    }
    get_subtree_size(root_node, subtree_size, n_children, children_of, n_edge);

    for (int32 node = n_tip; node != node_limit; node++) {
      quicksort_by_smallest(children_of[node],
                            children_of[node] + n_children[node] - 1,
                            subtree_size);
    }
    int32 * node_order = (int32*) malloc(n_node * sizeof(int32));
    for (int32 i = 0; i != n_node; ++i) {
      node_order[i] = i + n_tip;
    }
    quicksort_by_smallest(node_order, node_order + n_node - 1, subtree_size);
    std::free(subtree_size);

    IntegerMatrix ret(n_edge, 2);
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
      delete children_of[i];
    }
    delete[] (children_of);
    std::free(node_order);

    return (ret);
  }
}

#endif
