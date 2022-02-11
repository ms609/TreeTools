#include <Rcpp/Lightest>
#include "../inst/include/TreeTools.h"
using namespace Rcpp;

const uintx powers_of_two[16] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,
                                 2048, 4096, 8192, 16384, 32768};
const uintx BIN_SIZE = 8;

#define PO_PARENT(i) edge(order[i], 0)
#define PO_CHILD(i) edge(order[i], 1)


// Edges must be listed in 'strict' postorder, i.e. two-by-two
// [[Rcpp::export]]
RawMatrix cpp_edge_to_splits(const IntegerMatrix edge,
                             const IntegerVector order,
                             const IntegerVector nTip) {
  // Check input is valid
  if (edge.cols() != 2) {
    throw std::invalid_argument("Edge matrix must contain two columns");
  }
  if (1UL + edge.rows() > UINTX_MAX - 1U) { /* UINT_MAX denotes NOT_TRIVAL */
    throw std::length_error("Too many edges in tree for edge_to_splits: "       // # nocov
                            "Contact maintainer for advice");                   // # nocov
  }
  if (nTip[0] < 1) {
    throw(std::length_error("Tree must contain tips."));
  }

  // Initialize
  const uintx n_edge = edge.rows(),
              n_node = n_edge + 1,
              root_node = PO_PARENT(n_edge - 1),
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
  if (n_edge != order.length()) {
    stop("Length of `order` must equal number of edges");
  }

  uintx** splits = new uintx*[n_node];
  for (uintx i = 0; i != n_node; i++) {
    splits[i] = new uintx[n_bin](); // () zero-initializes
  }

  // Populate splits
  for (uintx i = 0; i != n_tip; ++i) {
    splits[i][uintx(i / BIN_SIZE)] = powers_of_two[i % BIN_SIZE];
  }

  uintx root_child = PO_CHILD(n_edge - 1);
  int32 root_children = 1;
  for (uintx i = 0; i != n_edge - 1; ++i) { // Omit last edge
    const uintx parent = PO_PARENT(i);
    const uintx child = PO_CHILD(i);
    if (parent == root_node) {
      ++root_children;
      if (child > n_tip) {
        root_child = uintx(child);
      }
    }
    for (uintx j = 0; j != n_bin; ++j) {
      splits[parent - 1][j] |= splits[child - 1][j];
    }
  }

  for (uintx i = 0; i != n_tip; ++i) {
    delete[] splits[i];
  }

  // Only return non-trivial splits
  uintx n_trivial = 0;
  const uintx NOT_TRIVIAL = UINTX_MAX,
              trivial_origin = root_node - 1,
              trivial_two = (root_children == 2 ? root_child - 1 : NOT_TRIVIAL),
              n_return = n_edge - n_tip - (trivial_two != NOT_TRIVIAL ? 1 : 0);
  RawMatrix ret(n_return, n_bin);
  IntegerVector names(n_return);

  for (uintx i = n_tip; i != n_node; ++i) {
    if (i == trivial_origin || i == trivial_two) {
      ++n_trivial;
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
