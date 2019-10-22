#include <Rcpp.h>
using namespace Rcpp;
#include <stdint.h> /* for uint32_t */

const uint32_t powers_of_two[32] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512,
                                    1024, 2048, 4096, 8192, 16384, 32768,
                                    2^16, 2^17, 2^18, 2^19, 2^20, 2^21, 2^22,
                                    2^23, 2^24, 2^25, 2^26, 2^27, 2^28, 2^29,
                                    2^30, 2^31};

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
      ret(i, j) = splits[n_tip + i + 1][j];
    }
  }

  return(ret);
}
