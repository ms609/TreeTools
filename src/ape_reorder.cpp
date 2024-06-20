#include <Rcpp/Lightest>
#include "ape_reorder.h"
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector ape_neworder_phylo(IntegerVector n_tips, IntegerVector parent,
                                 IntegerVector child, IntegerVector n_edges,
                                 IntegerVector order) {
  IntegerVector neworder(parent.size() * 2);
  
  ape_neworder_phylo(
    reinterpret_cast<int*>(n_tips.begin()),
    reinterpret_cast<int*>(parent.begin()),
    reinterpret_cast<int*>(child.begin()),
    reinterpret_cast<int*>(n_edges.begin()),
    reinterpret_cast<int*>(neworder.begin()),
    reinterpret_cast<int*>(order.begin())
  );
  
  return neworder;
}
  
// [[Rcpp::export]]
IntegerVector ape_neworder_pruningwise(
    IntegerVector n_tips, IntegerVector n_node, 
    IntegerVector parent, IntegerVector child,
    IntegerVector n_edges) {
  IntegerVector neworder(parent.size() * 2);
  
  ape_neworder_pruningwise(
    reinterpret_cast<int*>(n_tips.begin()),
    reinterpret_cast<int*>(n_node.begin()),
    reinterpret_cast<int*>(parent.begin()),
    reinterpret_cast<int*>(child.begin()),
    reinterpret_cast<int*>(n_edges.begin()),
    reinterpret_cast<int*>(neworder.begin())
  );
  
  return neworder;
}
