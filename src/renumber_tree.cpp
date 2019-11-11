#include <Rcpp.h>
#include <stdio.h>
using namespace Rcpp;

void report_calloc_error() {
  Rprintf("Error allocating memory in calloc");
}

/*
 * Returns the label of its parent, writing that parent's parent edge
 * if necessary
 */
void create_edge_leading_to(int old_child_label, int *parent_of,
                 int *new_label, int *next_label, int *next_edge,
                 int *final_parent, int *final_child) {

  if (!new_label[parent_of[old_child_label]]) {
    if (parent_of[parent_of[old_child_label]]) {
      /* Create the parent edge of the parent node */
      create_edge_leading_to(parent_of[old_child_label], parent_of,
                             new_label, next_label, next_edge,
                             final_parent, final_child);
    } else {
      /* Parent is the root node */
      new_label[parent_of[old_child_label]] = (*next_label)++;
    }
  }

  if (!new_label[old_child_label]) {
    new_label[old_child_label] = (*next_label)++;
  }
  final_parent[*next_edge] = new_label[parent_of[old_child_label]];
  final_child[*next_edge] = new_label[old_child_label];
  ++(*next_edge);

}

// [[Rcpp::export]]
IntegerMatrix order_edges_number_nodes(IntegerVector parent, IntegerVector child)
{
  const int n_edge = parent.length(), node_limit = n_edge + 2;
  if (child.length() != n_edge) {
    throw(std::length_error("Length of parent and child must match"));
  }

  int next_edge = 0,
    root_node = n_edge + n_edge;

  int * final_p =   (int*) calloc(n_edge, sizeof(int)),  /* calloc zero-initializes */
      * final_c =   (int*) calloc(n_edge, sizeof(int)),
      * parent_of = (int*) calloc(node_limit, sizeof(int)),
      * new_label = (int*) calloc(node_limit, sizeof(int));

  if (final_p == NULL) {
    report_calloc_error();
  } else if (final_c == NULL) {
    free(final_p);
    report_calloc_error();
  } else if (parent_of == NULL) {
    free(final_p);
    free(final_c);
    report_calloc_error();
  } else if (new_label == NULL) {
    free(final_p);
    free(final_c);
    free(parent_of);
    report_calloc_error();
  }

  for (int i = 0; i < n_edge; i++) {
    /* Initialize */
    parent_of[child[i]] = parent[i];
    if (parent[i] < root_node) root_node = parent[i];
  }
  parent_of[root_node] = 0;
  int next_label = root_node;

  for (int tip = 1; tip < root_node; tip++) {
    new_label[tip] = tip;
    create_edge_leading_to(tip, parent_of,
                           new_label, &next_label, &next_edge,
                           final_p, final_c);
  }

  IntegerMatrix ret(n_edge, 2);
  for (int i = 0; i < n_edge; i++) {
    ret(i, 0) = final_p[i];
    ret(i, 1) = final_c[i];
  }

  free(final_p);
  free(final_c);
  free(parent_of);
  free(new_label);
  return (ret);
}
