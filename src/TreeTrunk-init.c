#define USE_RINTERNALS

#include <Rinternals.h>
#include <stdlib.h> /* for NULL */

#include "ape_reorder.h"
#include "renumber_tree.h"

extern SEXP _TreeTrunk_phangorn_bipCPP(SEXP, SEXP);

static const R_CMethodDef cMethods[] = {
  {"order_edges_number_nodes", (DL_FUNC) &order_edges_number_nodes, 3, order_edges_number_nodes_t},
  {"ape_neworder_phylo",       (DL_FUNC) &ape_neworder_phylo, 6, ape_neworder_phylo_t},
  {"ape_node_depth",           (DL_FUNC) &ape_node_depth, 7, ape_node_depth_t},
  {"ape_neworder_pruningwise", (DL_FUNC) &ape_neworder_pruningwise, 6, ape_neworder_pruningwise_t},
  {NULL, NULL, 0, NULL}
};

static const R_CallMethodDef callMethods[] = {
  {"_TreeTrunk_phangorn_bipCPP",   (DL_FUNC) &_TreeTrunk_phangorn_bipCPP, 2},
  {"RENUMBER_TREE",  (DL_FUNC) &RENUMBER_TREE,  3},
  {"RENUMBER_EDGES", (DL_FUNC) &RENUMBER_EDGES, 3},
  {NULL, NULL, 0}
};

void R_init_TreeSearch(DllInfo *dll) {
  R_registerRoutines(dll, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
