#include <Rcpp.h>
#include "types.h"
using namespace Rcpp;

const uintx powers_of_two[16] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,
                                 2048, 4096, 8192, 16384, 32768};
const intx BIN_SIZE = 8;


// Edges must be listed in 'strict' postorder, i.e. two-by-two
// [[Rcpp::export]]
RawMatrix cpp_edge_to_splits(IntegerMatrix edge, IntegerVector nTip) {
  if (edge.cols() != 2) {
    throw std::invalid_argument("Edge matrix must contain two columns");
  }
  if (1L + edge.rows() > INTX_MAX) {
    throw(std::length_error("Too many edges in tree for edge_to_splits: "
                              "Contact maintainer for advice"));
  }

  const intx n_edge = edge.rows(),
             n_node = n_edge + 1,
             n_tip = nTip[0],
             n_bin = ((n_tip - 1) / BIN_SIZE) + 1;

  if (n_edge == n_tip) { /* No internal nodes resolved */
    return RawMatrix (0, n_bin);
  }
  if (n_edge < 3) {
    /* Cannot calculate trivial_two below. */
    throw(std::length_error("Not enough edges in tree for edge_to_splits."));
  }

  uintx** splits = new uintx*[n_node];
  for (intx i = 0; i != n_node; i++) {
    splits[i] = new uintx[n_bin];
    for (intx j = 0; j != n_bin; j++) {
      splits[i][j] = 0;
    }
  }

  for (intx i = 0; i != n_tip; i++) {
    splits[i][intx(i / BIN_SIZE)] = powers_of_two[i % BIN_SIZE];
  }

  for (intx i = 0; i != n_edge - 1; i++) { /* final edge is second root edge */
    for (intx j = 0; j != n_bin; j++) {
      splits[intx(edge(i, 0) - 1)][j] |= splits[intx(edge(i, 1) - 1)][j];
    }
  }

  for (intx i = 0; i != n_tip; i++) {
    delete[] splits[i];
  }
  
  intx n_trivial = 0;
  const intx NOT_TRIVIAL = -1;
  const intx trivial_origin = edge(n_edge - 1, 0) - 1,
    trivial_two = (edge(n_edge - 1, 0) == edge(n_edge - 3, 0) ?
                     NOT_TRIVIAL : (edge(n_edge - 1, 1) - 1L));
  const intx n_return = n_edge - n_tip - (trivial_two != NOT_TRIVIAL ? 1 : 0);
  RawMatrix ret(n_return, n_bin);
  IntegerVector names(n_return);

  for (intx i = n_tip; i != n_node; i++) {
    if (i == trivial_origin || i == trivial_two) {
      n_trivial++;
    } else {
      for (intx j = 0; j != n_bin; j++) {
        ret(i - n_tip - n_trivial, j) = splits[i][j];
        names[i - n_tip - n_trivial] = (i + 1);
      }
    }
    delete[] splits[i];
  }
  
  delete[] splits;

  rownames(ret) = names;
  return(ret);
}
