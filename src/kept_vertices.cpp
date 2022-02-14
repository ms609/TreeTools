#include <Rcpp/Lightest>
#include "../inst/include/TreeTools/types.h"
using namespace Rcpp;

// edge must be a two-column edge matrix in preorder
// [[Rcpp::export]]
IntegerVector kept_vertices (const IntegerMatrix edge,
                             const LogicalVector kept) {
  IntegerVector ret(edge.nrow() + 2);
  for (intx i = kept.length(); i--; ) {
    if (kept[i]) {
      ret[i + 1] = 2;
    }
  }
  for (intx i = edge.nrow(); i--; ) {
    if (ret[edge(i, 1)]) {
      ++ret[edge(i, 0)];
    }
  }
  return ret;
}