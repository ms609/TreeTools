#define USE_RINTERNALS

#include <Rinternals.h>
#include <stdlib.h> /* for NULL */

#include "ape_reorder.h"

extern SEXP _TreeTools_phangorn_bipCPP(SEXP, SEXP);
extern SEXP _TreeTools_cpp_edge_to_splits(SEXP, SEXP);
extern SEXP _TreeTools_preorder_edges_and_nodes(SEXP, SEXP);
extern SEXP _TreeTools_num_to_parent(SEXP, SEXP);
extern SEXP _TreeTools_edge_to_num(SEXP, SEXP, SEXP);
extern SEXP _TreeTools_edge_to_shape(SEXP, SEXP, SEXP);

static const R_CMethodDef cMethods[] = {
  {"ape_neworder_phylo",       (DL_FUNC) &ape_neworder_phylo, 6, ape_neworder_phylo_t},
  {"ape_node_depth",           (DL_FUNC) &ape_node_depth, 7, ape_node_depth_t},
  {"ape_neworder_pruningwise", (DL_FUNC) &ape_neworder_pruningwise, 6, ape_neworder_pruningwise_t},
  {NULL, NULL, 0, NULL}
};

static const R_CallMethodDef callMethods[] = {
  {"_TreeTools_num_to_parent", (DL_FUNC) &_TreeTools_num_to_parent, 2},
  {"_TreeTools_edge_to_num", (DL_FUNC) &_TreeTools_edge_to_num, 3},
  {"_TreeTools_edge_to_shape", (DL_FUNC) &_TreeTools_edge_to_shape, 3},
  {"_TreeTools_cpp_edge_to_splits", (DL_FUNC) &_TreeTools_cpp_edge_to_splits, 2},
  {"_TreeTools_phangorn_bipCPP", (DL_FUNC) &_TreeTools_phangorn_bipCPP, 2},
  {"_TreeTools_preorder_edges_and_nodes", (DL_FUNC) &_TreeTools_preorder_edges_and_nodes, 2},
  {NULL, NULL, 0}
};

void R_init_TreeTools(DllInfo *dll) {
  R_registerRoutines(dll, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
