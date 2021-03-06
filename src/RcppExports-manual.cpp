#include "../inst/include/TreeTools.h"
#include <Rcpp.h>

using namespace Rcpp;

// preorder_edges_and_nodes
RcppExport SEXP _TreeTools_preorder_edges_and_nodes(SEXP parentSEXP, SEXP childSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< const IntegerVector >::type parent(parentSEXP);
  Rcpp::traits::input_parameter< const IntegerVector >::type child(childSEXP);
  rcpp_result_gen = Rcpp::wrap(TreeTools::preorder_edges_and_nodes(parent, child));
  return rcpp_result_gen;
  END_RCPP
}
// postorder_edges
RcppExport SEXP _TreeTools_postorder_edges(SEXP edgeSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< const IntegerMatrix >::type edge(edgeSEXP);
  rcpp_result_gen = Rcpp::wrap(TreeTools::postorder_edges(edge));
  return rcpp_result_gen;
  END_RCPP
}
// root_binary
RcppExport SEXP _TreeTools_root_binary(SEXP edgeSEXP, SEXP outgroupSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< const IntegerMatrix >::type edge(edgeSEXP);
  Rcpp::traits::input_parameter< const int >::type outgroup(outgroupSEXP);
  rcpp_result_gen = Rcpp::wrap(TreeTools::root_binary(edge, outgroup));
  return rcpp_result_gen;
  END_RCPP
}
// root_on_node
RcppExport SEXP _TreeTools_root_on_node(SEXP phySEXP, SEXP outgroupSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< const List >::type phy(phySEXP);
  Rcpp::traits::input_parameter< const int >::type outgroup(outgroupSEXP);
  rcpp_result_gen = Rcpp::wrap(TreeTools::root_on_node(phy, outgroup));
  return rcpp_result_gen;
  END_RCPP
}
