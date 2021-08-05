#ifndef TreeTools_renumber_tree_
#define TreeTools_renumber_tree_

#include <cstdio>
#include <cstdlib> /* for calloc */
#include <Rcpp.h>
#include "types.h"
using namespace Rcpp;

namespace TreeTools {
  inline void swap(intx *a, intx *b) {
    const intx temp = *a;
    *a = *b;
    *b = temp;
  }

  /* Requires unsigned integers. */
  /* If we chose signed, we'd have to impose a limit on n_children, which
   * would exclude star trees */
  inline void quicksort_by_smallest(intx *left, const intx *right, const intx *sort_by) {
    if (left >= right) return;

    const intx pivot = sort_by[*right];
    intx *centre = left;
    for (intx *i = left; i <= right; i++) {
      if (sort_by[*i] <= pivot) {
        swap(centre, i);
        ++centre;
      }
    }
    quicksort_by_smallest(left, centre - 2, sort_by);
    quicksort_by_smallest(centre, right, sort_by);
  }

inline void add_child_edges(const intx node, const intx node_label,
                            intx const* const* children_of,
                            const intx *n_children,
                            IntegerMatrix final_edges,
                            intx *next_edge, intx *next_label,
                            const intx *n_tip, const intx *n_edge) {

    for (intx child = 0; child != n_children[node]; child++) {

      final_edges(*next_edge, 0) = node_label;
      intx this_child = children_of[node][child];

      if (this_child <= *n_tip) {

        final_edges(*next_edge, 1) = this_child;
        ++(*next_edge);

      } else {

        const intx child_label = (*next_label)++;
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
    if (2.0 * (2 + child.length()) > double(INTX_CONSERVATIVE_MAX)) {
      throw std::length_error("Too many edges in tree: "
                              "Contact 'TreeTools' maintainer for support.");
    }

    const intx n_edge = parent.length(),
               node_limit = n_edge + 2;

    if (child.length() != n_edge) {
      throw std::invalid_argument("Length of parent and child must match");
    }

    intx next_edge = 0,
         root_node = n_edge * 2, /* Initialize with too-big value */
         n_tip = 0;

    intx * parent_of = (intx*) std::calloc(node_limit, sizeof(intx)),
         * n_children = (intx*) std::calloc(node_limit, sizeof(intx)),
         * smallest_desc = (intx*) std::calloc(node_limit, sizeof(intx));
    intx ** children_of = new intx*[node_limit];

    for (intx i = n_edge; i--; ) {
      parent_of[child[i]] = parent[i];
      ++(n_children[parent[i]]);
    }

    for (intx i = 1; i != node_limit; i++) {
      if (!parent_of[i]) root_node = i;
      if (!n_children[i]) ++n_tip;
      children_of[i] = new intx[n_children[i]];
    }

    for (intx tip = 1; tip != n_tip + 1; ++tip) {
      smallest_desc[tip] = tip;
      intx parent = parent_of[tip];
      while (!smallest_desc[parent]) {
        smallest_desc[parent] = tip;
        parent = parent_of[parent];
      }
    }
    std::free(parent_of);

    intx * found_children = (intx*) std::calloc(node_limit, sizeof(intx));
    for (intx i = n_edge; i--; ) {
      children_of[parent[i]][(found_children[parent[i]])++] = child[i];
    }
    std::free(found_children);

    for (intx node = n_tip + 1; node != node_limit; node++) {
      quicksort_by_smallest(children_of[node],
                            children_of[node] + n_children[node] - 1,
                            smallest_desc);
    }
    std::free(smallest_desc);

    intx next_label = n_tip + 2;
    IntegerMatrix ret(n_edge, 2);
    add_child_edges(root_node, n_tip + 1,
                    children_of, n_children, ret,
                    &next_edge, &next_label, &n_tip, &n_edge);

    std::free(n_children);

    for (intx i = 1; i != node_limit; i++) {
      delete[] children_of[i];
    }
    delete[] children_of;

    return (ret);
  }

inline intx get_subtree_size(intx node, intx *subtree_size, intx *n_children,
                       intx **children_of, intx n_edge) {
    if (!subtree_size[node]) {
      for (intx i = n_children[node]; i--; ) {
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
    if (1.0 + edge.nrow() > double(INTX_CONSERVATIVE_MAX)) {
      throw std::length_error("Too many edges in tree for postorder_edges: "
                              "Contact maintainer for advice");
      // In theory we could use INTX_MAX, which is larger than 16 bits on linux,
      // or we could change intx to 32 bit.  The former has caused a seg fault
      // with invalid permissions on linux builds, possibly related to callocing
      // children_of?
      // Rather than attempt to debug now, I've chosen to place a hard limit on
      // edge.nrow for the time being.  --MS, 2020-05-26
      // Possibly fixed in branch optim-postorder, to check: 2021-08-05
    }

    const intx
      n_edge = edge.nrow(),
      node_limit = n_edge + 1;

    intx
      root_node = 0,
      n_tip = 0;

    // 0.9999 leaves room for memory overhead: seems in practice to avoid
    // attempting a doomed call to calloc.
    if (long(n_edge * node_limit * sizeof(intx)) > 0.9999L * INTPTR_MAX) {
      throw std::length_error("Tree too large for postorder_edges. "            // # nocov
                              "Try running 64-bit R?");                         // # nocov
    }

    intx * parent_of = (intx*) std::calloc(node_limit, sizeof(intx)),
         * n_children = (intx*) std::calloc(node_limit, sizeof(intx)),
         * subtree_size = (intx*) std::calloc(node_limit, sizeof(intx));
    intx ** children_of = new intx*[node_limit];
    
    for (intx i = n_edge; i--; ) {
      parent_of[edge(i, 1)] = edge(i, 0);
      ++(n_children[edge(i, 0)]);
    }

    for (intx i = node_limit; i--; ) {
      if (parent_of[i] == 0) root_node = i;
      if (n_children[i] == 0) ++n_tip;
      children_of[i] = new intx[n_children[i]];
    }
    std::free(parent_of);

    intx * found_children = (intx*) std::calloc(node_limit, sizeof(intx));
    for (intx i = n_edge; i--; ) {
      children_of[edge(i, 0)][(found_children[edge(i, 0)])++] = edge(i, 1);
    }
    std::free(found_children);
    
    const intx n_node = n_edge - n_tip + 1;

    for (intx tip = 0; tip != n_tip; tip++) {
      subtree_size[tip] = 1;
    }
    get_subtree_size(root_node, subtree_size, n_children, children_of, n_edge);

    for (intx node = n_tip; node != node_limit; node++) {
      quicksort_by_smallest(children_of[node],
                            children_of[node] + n_children[node] - 1,
                            subtree_size);
    }
    intx * node_order = (intx*) malloc(n_node * sizeof(intx));
    for (intx i = 0; i != n_node; ++i) {
      node_order[i] = i + n_tip;
    }
    quicksort_by_smallest(node_order, node_order + n_node - 1, subtree_size);
    std::free(subtree_size);

    IntegerMatrix ret(n_edge, 2);
    intx this_edge = 0;
    for (intx i = 0; i != n_node; ++i) {
      const intx this_parent = node_order[i];
      for (intx j = 0; j != n_children[this_parent]; ++j) {
        ret(this_edge, 0) = this_parent + 1;
        ret(this_edge, 1) = children_of[this_parent][j] + 1;
        ++this_edge;
      }
    }
    std::free(n_children);
    for (intx i = node_limit; i--; ) {
      delete children_of[i];
    }
    delete[] (children_of);
    std::free(node_order);

    return (ret);
  }
}

#endif
