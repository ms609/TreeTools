#include <Rcpp.h>
#include "../inst/include/TreeTools.h"
using namespace Rcpp;

const uintx powers_of_two[16] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,
                                 2048, 4096, 8192, 16384, 32768};
const uintx BIN_SIZE = 8;


// Edges must be listed in 'strict' postorder, i.e. two-by-two
// [[Rcpp::export]]
RawMatrix cpp_edge_to_splits(IntegerMatrix edge, IntegerVector nTip) {
  if (edge.cols() != 2) {
    throw std::invalid_argument("Edge matrix must contain two columns");
  }
  if (1UL + edge.rows() > UINTX_MAX - 1U) { /* UINT_MAX denotes NOT_TRIVAL */
    throw(std::length_error("Too many edges in tree for edge_to_splits: "       // # nocov
                            "Contact maintainer for advice"));                  // # nocov
  }
  if (nTip[0] < 1) {
    throw(std::length_error("Tree must contain tips."));
  }

  const uintx n_edge = edge.rows(),
              n_node = n_edge + 1,
              n_tip = nTip[0],
              n_bin = ((n_tip - 1) / BIN_SIZE) + 1;

  if (n_edge == n_tip  /* No internal nodes resolved */
        || n_tip < 4) { /* Need four tips to split non-trivially */
    return RawMatrix (0, n_bin);
  }
  if (n_edge < 3) {
    /* Cannot calculate trivial_two below. */
    throw(std::length_error("Not enough edges in tree for edge_to_splits."));
  }

  uintx** splits = new uintx*[n_node];
  for (uintx i = 0; i != n_node; i++) {
    splits[i] = new uintx[n_bin](); // () zero-initializes
  }

  for (uintx i = 0; i != n_tip; i++) {
    splits[i][uintx(i / BIN_SIZE)] = powers_of_two[i % BIN_SIZE];
  }

  for (uintx i = 0; i != n_edge - 1; i++) { /* final edge is second root edge */
    for (uintx j = 0; j != n_bin; j++) {
      splits[uintx(edge(i, 0) - 1)][j] |= splits[uintx(edge(i, 1) - 1)][j];
    }
  }

  for (uintx i = 0; i != n_tip; i++) {
    delete[] splits[i];
  }

  uintx n_trivial = 0;
  const uintx NOT_TRIVIAL = UINTX_MAX,
              trivial_origin = edge(n_edge - 1, 0) - 1,
              trivial_two = (edge(n_edge - 1, 0) == edge(n_edge - 3, 0) ?
                            NOT_TRIVIAL : (edge(n_edge - 1, 1) - 1L)),
              n_return = n_edge - n_tip - (trivial_two != NOT_TRIVIAL ? 1 : 0);
  RawMatrix ret(n_return, n_bin);
  IntegerVector names(n_return);

  for (uintx i = n_tip; i != n_node; i++) {
    if (i == trivial_origin || i == trivial_two) {
      n_trivial++;
    } else {
      for (uintx j = 0; j != n_bin; j++) {
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
