#include <cstdio>
#include <cstdlib> /* for calloc */
#include <Rcpp.h>
using namespace Rcpp;

int smallest_descendant(int *node, int *smallest_desc, int *n_children,
                        int *children_of, const int *n_edge) {
  if (!smallest_desc[*node]) {
    smallest_desc[*node] =
      smallest_descendant(&children_of[*node * *n_edge + 0], smallest_desc,
                          n_children, children_of, n_edge);
    for (int j = 1; j < n_children[*node]; j++) {
      int this_child =
        smallest_descendant(&children_of[*node * *n_edge + j], smallest_desc,
                            n_children, children_of, n_edge);
      smallest_desc[*node] = smallest_desc[*node] < this_child
        ? smallest_desc[*node] : this_child;
    }
  }
  return smallest_desc[*node];
}

void swap(int *a, int *b) {
  const int temp = *a;
  *a = *b;
  *b = temp;
}

void quicksort_by_smallest(int *to_sort, int *sort_by, int left, int right) {
  if (left >= right) return;
  int pivot = sort_by[to_sort[right]];
  int centre = left;
  for (int i = left; i <= right; i++) {
    if (sort_by[to_sort[i]] <= pivot) {
      swap(&to_sort[centre], &to_sort[i]);
      centre++;
    }
  }
  quicksort_by_smallest(to_sort, sort_by, left, centre - 2);
  quicksort_by_smallest(to_sort, sort_by, centre, right);
}

void add_child_edges(int node, int node_label,
                     int *children_of, int *n_children,
                     IntegerMatrix final_edges, int *next_edge, int *next_label,
                     int *n_tip, const int *n_edge) {

  for (int child = 0; child < n_children[node]; child++) {

    final_edges(*next_edge, 0) = node_label;
    int this_child = children_of[node * *n_edge + child];

    if (this_child <= *n_tip) {

      final_edges(*next_edge, 1) = this_child;
      ++(*next_edge);

    } else {

      int child_label = (*next_label)++;
      final_edges(*next_edge, 1) = child_label;
      ++(*next_edge);

      add_child_edges(this_child, child_label, children_of, n_children,
                      final_edges, next_edge, next_label, n_tip, n_edge);

    }
  }
}

// [[Rcpp::export]]
IntegerMatrix preorder_edges_and_nodes(IntegerVector parent, IntegerVector child)
{
  const int n_edge = parent.length(), node_limit = n_edge + 2;
  if (child.length() != n_edge) {
    throw(std::length_error("Length of parent and child must match"));
  }

  int next_edge = 0,
    root_node = n_edge + n_edge,
    n_tip = 0;

  int * parent_of = (int*) std::calloc(node_limit, sizeof(int)),
      * n_children = (int*) std::calloc(node_limit, sizeof(int)),
      * smallest_desc = (int*) std::calloc(node_limit, sizeof(int)),
      * children_of = (int*) std::calloc(n_edge * node_limit, sizeof(int));

  for (int i = 0; i < n_edge; i++) {
    parent_of[child[i]] = parent[i];
    children_of[parent[i] * n_edge + n_children[parent[i]]] = child[i];
    (n_children[parent[i]])++;
  }

  for (int i = 1; i < node_limit; i++) {
    if (parent_of[i] == 0) root_node = i;
    if (n_children[i] == 0) ++n_tip;
  }
  std::free(parent_of);

  for (int tip = 1; tip <= n_tip; tip++) {
    smallest_desc[tip] = tip;
  }

  for (int node = n_tip + 1; node < node_limit; node++) {
    smallest_descendant(&node, smallest_desc, n_children, children_of, &n_edge);
    quicksort_by_smallest(&children_of[node * n_edge], smallest_desc,
                          0, n_children[node] - 1);
  }
  std::free(smallest_desc);

  int next_label = n_tip + 2;
  IntegerMatrix ret(n_edge, 2);
  add_child_edges(root_node, n_tip + 1,
                  children_of, n_children, ret,
                  &next_edge, &next_label, &n_tip, &n_edge);

  std::free(n_children);
  std::free(children_of);

  return (ret);
}

int get_subtree_size(int node, int *subtree_size, int *n_children,
                     int *children_of, int n_edge) {
  if (!subtree_size[node]) {
    for (int i = 0; i < n_children[node]; i++) {
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
  const int n_edge = edge.nrow(), node_limit = n_edge + 1;
  int root_node = 0, n_tip = 0;
  int * parent_of = (int*) std::calloc(node_limit, sizeof(int)),
      * n_children = (int*) std::calloc(node_limit, sizeof(int)),
      * subtree_size = (int*) std::calloc(node_limit, sizeof(int)),
      * children_of = (int*) std::calloc(n_edge * node_limit, sizeof(int));

  for (int i = 0; i < n_edge; i++) {
    parent_of[edge(i, 1)] = edge(i, 0);
    children_of[edge(i, 0) * n_edge + n_children[edge(i, 0)]] = edge(i, 1);
    (n_children[edge(i, 0)])++;
  }

  for (int i = 0; i < node_limit; i++) {
    if (parent_of[i] == 0) root_node = i;
    if (n_children[i] == 0) ++n_tip;
  }
  std::free(parent_of);
  const int n_node = n_edge - n_tip + 1;

  for (int tip = 0; tip < n_tip; tip++) {
    subtree_size[tip] = 1;
  }
  get_subtree_size(root_node, subtree_size, n_children, children_of, n_edge);

  for (int node = n_tip; node < node_limit; node++) {
    quicksort_by_smallest(&children_of[node * n_edge], subtree_size,
                          0, n_children[node] - 1);
  }
  int * node_order = (int*) malloc(n_node * sizeof(int));
  for (int i = 0; i < n_node; i++) {
    node_order[i] = i + n_tip;
  }
  quicksort_by_smallest(node_order, subtree_size, 0, n_node - 1);

  IntegerMatrix ret(n_edge, 2);
  int this_edge = 0;
  for (int i = 0; i < n_node; i++) {
    const int this_parent = node_order[i];
    for (int j = 0; j < n_children[this_parent]; j++) {
      ret(this_edge, 0) = this_parent + 1;
      ret(this_edge, 1) = children_of[this_parent * n_edge + j] + 1;
      ++this_edge;
    }
  }
  std::free(n_children);
  std::free(subtree_size);
  std::free(children_of);
  std::free(node_order);

  return (ret);
}
