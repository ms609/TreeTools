#define USE_RINTERNALS

#include <Rinternals.h>
#include <stdlib.h> /* for NULL */

#include "ape_reorder.h"

extern SEXP _TreeTools_cpp_edge_to_splits(SEXP, SEXP);
extern SEXP _TreeTools_preorder_edges_and_nodes(SEXP, SEXP);
extern SEXP _TreeTools_postorder_edges(SEXP);
extern SEXP _TreeTools_root_on_node(SEXP, SEXP);
extern SEXP _TreeTools_root_binary(SEXP, SEXP);
extern SEXP _TreeTools_num_to_parent(SEXP, SEXP);
extern SEXP _TreeTools_edge_to_num(SEXP, SEXP, SEXP);
extern SEXP _TreeTools_edge_to_rooted_shape(SEXP, SEXP, SEXP);
extern SEXP _TreeTools_random_parent(SEXP, SEXP);
extern SEXP _TreeTools_rooted_shape_to_edge(SEXP, SEXP);
extern SEXP _TreeTools_tips_in_splits(SEXP, SEXP);
extern SEXP _TreeTools_as_newick(SEXP);

static const R_CMethodDef cMethods[] = {
  {"ape_neworder_phylo",       (DL_FUNC) &ape_neworder_phylo, 6, ape_neworder_phylo_t},
  {"ape_node_depth",           (DL_FUNC) &ape_node_depth, 7, ape_node_depth_t},
  {"ape_neworder_pruningwise", (DL_FUNC) &ape_neworder_pruningwise, 6, ape_neworder_pruningwise_t},
  {NULL, NULL, 0, NULL}
};

static const R_CallMethodDef callMethods[] = {
  {"_TreeTools_num_to_parent", (DL_FUNC) &_TreeTools_num_to_parent, 2},
  {"_TreeTools_edge_to_num", (DL_FUNC) &_TreeTools_edge_to_num, 3},
  {"_TreeTools_edge_to_rooted_shape", (DL_FUNC) &_TreeTools_edge_to_rooted_shape, 3},
  {"_TreeTools_random_parent", (DL_FUNC) &_TreeTools_random_parent, 2},
  {"_TreeTools_rooted_shape_to_edge", (DL_FUNC) &_TreeTools_rooted_shape_to_edge, 2},
  {"_TreeTools_cpp_edge_to_splits", (DL_FUNC) &_TreeTools_cpp_edge_to_splits, 2},
  {"_TreeTools_preorder_edges_and_nodes", (DL_FUNC) &_TreeTools_preorder_edges_and_nodes, 2},
  {"_TreeTools_postorder_edges", (DL_FUNC) &_TreeTools_postorder_edges, 1},
  {"_TreeTools_root_on_node", (DL_FUNC) &_TreeTools_root_on_node, 2},
  {"_TreeTools_root_binary", (DL_FUNC) &_TreeTools_root_binary, 2},
  {"_TreeTools_tips_in_splits", (DL_FUNC) &_TreeTools_tips_in_splits, 1},
  {"_TreeTools_as_newick", (DL_FUNC) &_TreeTools_as_newick, 1},
  {NULL, NULL, 0}
};

void R_init_TreeTools(DllInfo *dll) {
  R_registerRoutines(dll, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
