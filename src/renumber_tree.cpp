#include <cstdio>
#include <cstdlib> /* for calloc */
#include <Rcpp.h>
#include "types.h"
using namespace Rcpp;

xint smallest_descendant(xint *node, xint *smallest_desc, xint *n_children,
                         xint *children_of, const xint *n_edge) {
  if (!smallest_desc[*node]) {
    smallest_desc[*node] =
      smallest_descendant(&children_of[*node * *n_edge + 0], smallest_desc,
                          n_children, children_of, n_edge);
    for (xint j = 1; j < n_children[*node]; j++) {
      xint this_child =
        smallest_descendant(&children_of[*node * *n_edge + j], smallest_desc,
                            n_children, children_of, n_edge);
      smallest_desc[*node] = smallest_desc[*node] < this_child
        ? smallest_desc[*node] : this_child;
    }
  }
  return smallest_desc[*node];
}

void swap(xint *a, xint *b) {
  const xint temp = *a;
  *a = *b;
  *b = temp;
}

void quicksort_by_smallest(xint *to_sort, xint *sort_by, xint left, xint right) {
  if (left >= right) return;
  xint pivot = sort_by[to_sort[right]];
  xint centre = left;
  for (xint i = left; i <= right; i++) {
    if (sort_by[to_sort[i]] <= pivot) {
      swap(&to_sort[centre], &to_sort[i]);
      centre++;
    }
  }
  quicksort_by_smallest(to_sort, sort_by, left, centre - 2);
  quicksort_by_smallest(to_sort, sort_by, centre, right);
}

void add_child_edges(xint node, xint node_label,
                     xint *children_of, xint *n_children,
                     IntegerMatrix final_edges, xint *next_edge, xint *next_label,
                     xint *n_tip, const xint *n_edge) {

  for (xint child = 0; child < n_children[node]; child++) {

    final_edges(*next_edge, 0) = node_label;
    xint this_child = children_of[node * *n_edge + child];

    if (this_child <= *n_tip) {

      final_edges(*next_edge, 1) = this_child;
      ++(*next_edge);

    } else {

      xint child_label = (*next_label)++;
      final_edges(*next_edge, 1) = child_label;
      ++(*next_edge);

      add_child_edges(this_child, child_label, children_of, n_children,
                      final_edges, next_edge, next_label, n_tip, n_edge);

    }
  }
}

// [[Rcpp::export]]
IntegerMatrix preorder_edges_and_nodes(IntegerVector parent,
                                       IntegerVector child)
{
  const xint n_edge = parent.length(),
             node_limit = n_edge + 2;

  if (2L * (2L + child.length()) > XINT_MAX) {
    throw(std::length_error("Too many edges: Contact maintainer for support."));
  }
  if (child.length() != n_edge) {
    throw(std::invalid_argument("Length of parent and child must match"));
  }

  xint next_edge = 0,
       root_node = n_edge + n_edge, /* Initialize with too-big value */
       n_tip = 0;

  xint * parent_of = (xint*) std::calloc(node_limit, sizeof(xint)),
       * n_children = (xint*) std::calloc(node_limit, sizeof(xint)),
       * smallest_desc = (xint*) std::calloc(node_limit, sizeof(xint)),
       * children_of = (xint*) std::calloc(n_edge * node_limit, sizeof(xint));

  for (xint i = 0; i < n_edge; i++) {
    parent_of[child[i]] = parent[i];
    children_of[parent[i] * n_edge + n_children[parent[i]]] = child[i];
    (n_children[parent[i]])++;
  }

  for (xint i = 1; i < node_limit; i++) {
    if (parent_of[i] == 0) root_node = i;
    if (n_children[i] == 0) ++n_tip;
  }
  std::free(parent_of);

  for (xint tip = 1; tip <= n_tip; tip++) {
    smallest_desc[tip] = tip;
  }

  for (xint node = n_tip + 1; node < node_limit; node++) {
    smallest_descendant(&node, smallest_desc, n_children, children_of, &n_edge);
    quicksort_by_smallest(&children_of[node * n_edge], smallest_desc,
                          0, n_children[node] - 1);
  }
  std::free(smallest_desc);

  xint next_label = n_tip + 2;
  IntegerMatrix ret(n_edge, 2);
  add_child_edges(root_node, n_tip + 1,
                  children_of, n_children, ret,
                  &next_edge, &next_label, &n_tip, &n_edge);

  std::free(n_children);
  std::free(children_of);

  return (ret);
}

xint get_subtree_size(xint node, xint *subtree_size, xint *n_children,
                     xint *children_of, xint n_edge) {
  if (!subtree_size[node]) {
    for (xint i = 0; i < n_children[node]; i++) {
      subtree_size[node] += get_subtree_size(children_of[node * n_edge + i],
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
IntegerMatrix postorder_edges(IntegerMatrix edge)
{
  if (1L + edge.nrow() > XINT_MAX) {
    throw(std::length_error("Too many edges in tree for postorder_edges: "
                            "Contact maintainer for advice"));
  }
  const xint n_edge = edge.nrow(), node_limit = n_edge + 1;
  xint root_node = 0, n_tip = 0;
  xint * parent_of = (xint*) std::calloc(node_limit, sizeof(xint)),
       * n_children = (xint*) std::calloc(node_limit, sizeof(xint)),
       * subtree_size = (xint*) std::calloc(node_limit, sizeof(xint)),
       * children_of = (xint*) std::calloc(n_edge * node_limit, sizeof(xint));

  for (xint i = 0; i < n_edge; i++) {
    parent_of[edge(i, 1)] = edge(i, 0);
    children_of[edge(i, 0) * n_edge + n_children[edge(i, 0)]] = edge(i, 1);
    (n_children[edge(i, 0)])++;
  }

  for (xint i = 0; i < node_limit; i++) {
    if (parent_of[i] == 0) root_node = i;
    if (n_children[i] == 0) ++n_tip;
  }
  std::free(parent_of);

  const xint n_node = n_edge - n_tip + 1;

  for (xint tip = 0; tip < n_tip; tip++) {
    subtree_size[tip] = 1;
  }
  get_subtree_size(root_node, subtree_size, n_children, children_of, n_edge);

  for (xint node = n_tip; node < node_limit; node++) {
    quicksort_by_smallest(&children_of[node * n_edge], subtree_size,
                          0, n_children[node] - 1);
  }
  xint * node_order = (xint*) malloc(n_node * sizeof(xint));
  for (xint i = 0; i < n_node; i++) {
    node_order[i] = i + n_tip;
  }
  quicksort_by_smallest(node_order, subtree_size, 0, n_node - 1);
  std::free(subtree_size);

  IntegerMatrix ret(n_edge, 2);
  xint this_edge = 0;
  for (xint i = 0; i < n_node; i++) {
    const xint this_parent = node_order[i];
    for (xint j = 0; j < n_children[this_parent]; j++) {
      ret(this_edge, 0) = this_parent + 1;
      ret(this_edge, 1) = children_of[this_parent * n_edge + j] + 1;
      ++this_edge;
    }
  }
  std::free(n_children);
  std::free(children_of);
  std::free(node_order);

  return (ret);
}
