#include <Rcpp.h>
#include <stdio.h>
using namespace Rcpp;

void report_calloc_error() {
  Rprintf("Error allocating memory in calloc");
}

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

void quicksort_by_smallest(int *to_sort, int *smallest_desc, int left, int right) {
  if (left >= right) return;
  int pivot = smallest_desc[to_sort[right]];
  int centre = left;
  for (int i = left; i <= right; i++) {
    if (smallest_desc[to_sort[i]] <= pivot) {
      swap(&to_sort[centre], &to_sort[i]);
      centre++;
    }
  }
  quicksort_by_smallest(to_sort, smallest_desc, left, centre - 2);
  quicksort_by_smallest(to_sort, smallest_desc, centre, right);
}

void add_child_edges(int node, int node_label,
                     int *children_of, int *n_children,
                     int *final_parent, int *final_child,
                     int *next_edge, int *next_label,
                     int *n_tip, const int *n_edge) {
  for (int child = 0; child < n_children[node]; child++) {
    final_parent[*next_edge] = node_label;
    int this_child = children_of[node * *n_edge + child];
    if (this_child <= *n_tip) {
      final_child[*next_edge] = this_child;
      ++(*next_edge);
    } else {
      int child_label = (*next_label)++;
      final_child[*next_edge] = child_label;
      ++(*next_edge);
      add_child_edges(this_child, child_label, children_of, n_children,
                      final_parent, final_child,
                      next_edge, next_label,
                      n_tip, n_edge);
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

  int * final_p =   (int*) calloc(n_edge, sizeof(int)),  /* calloc zero-initializes */
      * final_c =   (int*) calloc(n_edge, sizeof(int)),
      * parent_of = (int*) calloc(node_limit, sizeof(int)),
      * n_children = (int*) calloc(node_limit, sizeof(int)),
      * smallest_desc = (int*) calloc(node_limit, sizeof(int)),
      * children_of = (int*) calloc(n_edge * node_limit, sizeof(int));

  if (final_p == NULL || final_c == NULL || parent_of == NULL
        || n_children == NULL || smallest_desc == NULL || children_of == NULL) {
    /* Should probably free those that have been set */
    report_calloc_error();
  } /* Assume that calloc's fine otherwise... lazy but simple */


  for (int i = 0; i < n_edge; i++) {
    parent_of[child[i]] = parent[i];
    children_of[parent[i] * n_edge + n_children[parent[i]]] = child[i];
    (n_children[parent[i]])++;
  }

  for (int i = 1; i < node_limit; i++) {
    if (parent_of[i] == 0) root_node = i;
    if (n_children[i] == 0) ++n_tip;
  }

  for (int tip = 1; tip <= n_tip; tip++) {
    smallest_desc[tip] = tip;
  }

  for (int node = n_tip + 1; node < node_limit; node++) {
    smallest_descendant(&node, smallest_desc, n_children, children_of, &n_edge);
    quicksort_by_smallest(&children_of[node * n_edge], smallest_desc,
                          0, n_children[node] - 1);
  }

  int next_label = n_tip + 2;
  add_child_edges(root_node, n_tip + 1,
                  children_of, n_children,
                  final_p, final_c,
                  &next_edge, &next_label, &n_tip, &n_edge);

  IntegerMatrix ret(n_edge, 2);
  for (int i = 0; i < n_edge; i++) {
    ret(i, 0) = final_p[i];
    ret(i, 1) = final_c[i];
  }

  free(final_p);
  free(final_c);
  free(parent_of);
  free(n_children);
  free(smallest_desc);
  free(children_of);
  return (ret);
}
