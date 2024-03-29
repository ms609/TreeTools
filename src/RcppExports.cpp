// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/TreeTools.h"
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ClusterTable_new
SEXP ClusterTable_new(Rcpp::List phylo);
RcppExport SEXP _TreeTools_ClusterTable_new(SEXP phyloSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type phylo(phyloSEXP);
    rcpp_result_gen = Rcpp::wrap(ClusterTable_new(phylo));
    return rcpp_result_gen;
END_RCPP
}
// ClusterTable_matrix
Rcpp::IntegerMatrix ClusterTable_matrix(SEXP xp);
RcppExport SEXP _TreeTools_ClusterTable_matrix(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(ClusterTable_matrix(xp));
    return rcpp_result_gen;
END_RCPP
}
// ClusterTable_decode
Rcpp::IntegerVector ClusterTable_decode(SEXP xp);
RcppExport SEXP _TreeTools_ClusterTable_decode(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(ClusterTable_decode(xp));
    return rcpp_result_gen;
END_RCPP
}
// as_newick
CharacterVector as_newick(const IntegerMatrix edge);
RcppExport SEXP _TreeTools_as_newick(SEXP edgeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type edge(edgeSEXP);
    rcpp_result_gen = Rcpp::wrap(as_newick(edge));
    return rcpp_result_gen;
END_RCPP
}
// consensus_tree
LogicalMatrix consensus_tree(const List trees, const NumericVector p);
RcppExport SEXP _TreeTools_consensus_tree(SEXP treesSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type trees(treesSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(consensus_tree(trees, p));
    return rcpp_result_gen;
END_RCPP
}
// num_to_parent
IntegerVector num_to_parent(const IntegerVector n, const IntegerVector nTip);
RcppExport SEXP _TreeTools_num_to_parent(SEXP nSEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(num_to_parent(n, nTip));
    return rcpp_result_gen;
END_RCPP
}
// random_parent
IntegerVector random_parent(const IntegerVector nTip, const IntegerVector seed);
RcppExport SEXP _TreeTools_random_parent(SEXP nTipSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type nTip(nTipSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(random_parent(nTip, seed));
    return rcpp_result_gen;
END_RCPP
}
// edge_to_num
IntegerVector edge_to_num(IntegerVector parent, IntegerVector child, IntegerVector nTip);
RcppExport SEXP _TreeTools_edge_to_num(SEXP parentSEXP, SEXP childSEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type parent(parentSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type child(childSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(edge_to_num(parent, child, nTip));
    return rcpp_result_gen;
END_RCPP
}
// edge_to_mixed_base
IntegerVector edge_to_mixed_base(const IntegerVector parent, const IntegerVector child, const IntegerVector nTip);
RcppExport SEXP _TreeTools_edge_to_mixed_base(SEXP parentSEXP, SEXP childSEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type parent(parentSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type child(childSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(edge_to_mixed_base(parent, child, nTip));
    return rcpp_result_gen;
END_RCPP
}
// mixed_base_to_parent
IntegerVector mixed_base_to_parent(const IntegerVector n, const IntegerVector nTip);
RcppExport SEXP _TreeTools_mixed_base_to_parent(SEXP nSEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(mixed_base_to_parent(n, nTip));
    return rcpp_result_gen;
END_RCPP
}
// kept_vertices
IntegerVector kept_vertices(const IntegerMatrix edge, const LogicalVector kept);
RcppExport SEXP _TreeTools_kept_vertices(SEXP edgeSEXP, SEXP keptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type edge(edgeSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type kept(keptSEXP);
    rcpp_result_gen = Rcpp::wrap(kept_vertices(edge, kept));
    return rcpp_result_gen;
END_RCPP
}
// minimum_spanning_tree
IntegerMatrix minimum_spanning_tree(const IntegerVector order);
RcppExport SEXP _TreeTools_minimum_spanning_tree(SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(minimum_spanning_tree(order));
    return rcpp_result_gen;
END_RCPP
}
// path_lengths
NumericMatrix path_lengths(const IntegerMatrix edge, const DoubleVector weight);
RcppExport SEXP _TreeTools_path_lengths(SEXP edgeSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type edge(edgeSEXP);
    Rcpp::traits::input_parameter< const DoubleVector >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(path_lengths(edge, weight));
    return rcpp_result_gen;
END_RCPP
}
// cpp_edge_to_splits
RawMatrix cpp_edge_to_splits(const IntegerMatrix edge, const IntegerVector order, const IntegerVector nTip);
RcppExport SEXP _TreeTools_cpp_edge_to_splits(SEXP edgeSEXP, SEXP orderSEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type edge(edgeSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type order(orderSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_edge_to_splits(edge, order, nTip));
    return rcpp_result_gen;
END_RCPP
}
// duplicated_splits
LogicalVector duplicated_splits(const RawMatrix splits, const LogicalVector fromLast);
RcppExport SEXP _TreeTools_duplicated_splits(SEXP splitsSEXP, SEXP fromLastSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawMatrix >::type splits(splitsSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type fromLast(fromLastSEXP);
    rcpp_result_gen = Rcpp::wrap(duplicated_splits(splits, fromLast));
    return rcpp_result_gen;
END_RCPP
}
// mask_splits
RawMatrix mask_splits(RawMatrix x);
RcppExport SEXP _TreeTools_mask_splits(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RawMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(mask_splits(x));
    return rcpp_result_gen;
END_RCPP
}
// not_splits
RawMatrix not_splits(const RawMatrix x);
RcppExport SEXP _TreeTools_not_splits(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(not_splits(x));
    return rcpp_result_gen;
END_RCPP
}
// xor_splits
RawMatrix xor_splits(const RawMatrix x, const RawMatrix y);
RcppExport SEXP _TreeTools_xor_splits(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< const RawMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(xor_splits(x, y));
    return rcpp_result_gen;
END_RCPP
}
// and_splits
RawMatrix and_splits(const RawMatrix x, const RawMatrix y);
RcppExport SEXP _TreeTools_and_splits(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< const RawMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(and_splits(x, y));
    return rcpp_result_gen;
END_RCPP
}
// or_splits
RawMatrix or_splits(const RawMatrix x, const RawMatrix y);
RcppExport SEXP _TreeTools_or_splits(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< const RawMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(or_splits(x, y));
    return rcpp_result_gen;
END_RCPP
}
// thin_splits
RawMatrix thin_splits(const RawMatrix splits, const LogicalVector drop);
RcppExport SEXP _TreeTools_thin_splits(SEXP splitsSEXP, SEXP dropSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawMatrix >::type splits(splitsSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type drop(dropSEXP);
    rcpp_result_gen = Rcpp::wrap(thin_splits(splits, drop));
    return rcpp_result_gen;
END_RCPP
}
// splits_to_edge
IntegerMatrix splits_to_edge(const RawMatrix splits, const IntegerVector nTip);
RcppExport SEXP _TreeTools_splits_to_edge(SEXP splitsSEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawMatrix >::type splits(splitsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(splits_to_edge(splits, nTip));
    return rcpp_result_gen;
END_RCPP
}
// tips_in_splits
IntegerVector tips_in_splits(RawMatrix splits);
RcppExport SEXP _TreeTools_tips_in_splits(SEXP splitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RawMatrix >::type splits(splitsSEXP);
    rcpp_result_gen = Rcpp::wrap(tips_in_splits(splits));
    return rcpp_result_gen;
END_RCPP
}
// edge_to_rooted_shape
IntegerVector edge_to_rooted_shape(IntegerVector parent, IntegerVector child, IntegerVector nTip);
RcppExport SEXP _TreeTools_edge_to_rooted_shape(SEXP parentSEXP, SEXP childSEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type parent(parentSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type child(childSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(edge_to_rooted_shape(parent, child, nTip));
    return rcpp_result_gen;
END_RCPP
}
// rooted_shape_to_edge
IntegerMatrix rooted_shape_to_edge(NumericVector shape, IntegerVector nTip);
RcppExport SEXP _TreeTools_rooted_shape_to_edge(SEXP shapeSEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(rooted_shape_to_edge(shape, nTip));
    return rcpp_result_gen;
END_RCPP
}
