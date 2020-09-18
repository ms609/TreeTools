#ifndef TreeTools_H_GEN_
#define TreeTools_H_GEN_

#include <Rcpp.h>
using namespace Rcpp;

namespace TreeTools {

  using namespace Rcpp;

  extern IntegerMatrix preorder_edges_and_nodes(const IntegerVector,
                                                const IntegerVector);
  extern IntegerMatrix postorder_edges(const IntegerMatrix);
  extern IntegerMatrix root_on_node(const IntegerMatrix, const int);
}

#endif // TreeTools_H_GEN_
