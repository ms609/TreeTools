#include <Rcpp.h>
using namespace Rcpp;
#include <stdint.h> /* for uint32_t */

const uint32_t powers_of_two[32] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512,
                                    1024, 2048, 4096, 8192, 16384, 32768,
                                    65536, 131072, 262144, 524288, 1048576,
                                    2097152, 4194304, 8388608, 16777216,
                                    33554432, 67108864, 134217728, 268435456,
                                    536870912, 1073741824, 2147483648};

// [[Rcpp::export]]
NumericMatrix cpp_edge_to_splits(NumericMatrix edge) {
  if (edge.cols() != 2) {
    throw std::invalid_argument("Edge matrix must contain two columns");
  }

  const int n_edge = edge.rows(),
    n_node = n_edge + 1,
    n_tip = edge(0, 0) - 1,
    n_bin = n_tip / 32 + 1;

  uint32_t** splits = new uint32_t*[n_node];
  for (int i = 0; i < n_node; i++) {
    splits[i] = new uint32_t[n_bin];
    for (int j = 0; j < n_bin; j++) {
      splits[i][j] = 0;
    }
  }

  for (int i = 0; i < n_tip; i++) {
    splits[i][(int) i / 32] = powers_of_two[i % 32];
  }

  for (int i = n_edge - 1; i > 0; i--) { /* edge 0 is second root edge */
    int child = edge(i, 1) - 1;
    Rcout << "Edge " << i << ": parent = " << (edge(i, 0) - 1) << "."
          << splits[(int) (edge(i, 0) - 1)][0]
          << "; child = " << child
          << "." << splits[child][0] << "\n";
    for (int j = 0; j < n_bin; j++) {
      Rcout << "  Edge " << i << ", bin " << j
            << ": " << splits[(int) edge(i, 0) - 1][j] << " |= "
            << splits[child][j] << " = "
            << (splits[(int) edge(i, 0) - 1][j] | splits[child][j]) << ".\n";
      splits[(int) edge(i, 0) - 1][j] |= splits[child][j];
    }
  }

  NumericMatrix ret(n_edge - n_tip, n_bin);
  for (int i = 0; i < n_edge - n_tip; i++) {
    for (int j = 0; j < n_bin; j++) {
      ret(i, j) = splits[n_tip + i + 1][n_bin - j - 1];
    }
  }

  return(ret);
}
