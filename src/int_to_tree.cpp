#include <stdint.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector int_to_parent(NumericVector n, IntegerVector nTip) {
  if (nTip[0] < 2) {
    throw std::range_error("nTip must be > 1");
  }
  const unsigned int n_tip = nTip[0],
                     root_node = n_tip + n_tip - 1,
                     c_to_r = 1,
                     prime = n_tip - 2;
  uint64_t tree_id = n[0];
  unsigned int base;

  IntegerVector edge(n_tip + n_tip - 2);
  edge(0) = root_node;
  edge(1) = root_node;

  for (unsigned int i = 2; i < n_tip; i++) {
    base = (i + i - 3);
    const int i_prime = i + prime,
      i_prime_r = i_prime + c_to_r;

    uint64_t where = (tree_id % base) + 1;
    if (where >= i) {
      where += prime + 2 - i;
    }
    /*Rcout << i << ": " << tree_id << "%" << base << " => " << where << "\n";*/

    /*Rcout  << "   Re-parenting edge " << i_prime << " to parent of " << where
             << ", " << edge(where) << "\n";*/
    edge(i_prime) = edge(where);
    /*Rcout  << "   Re-parenting edge " << i << " to " << i_prime_r << "\n";*/
    edge(i) = i_prime_r;
    /*Rcout  << "   Re-parenting edge " << where << " to " << i_prime_r << "\n";*/
    edge(where) = i_prime_r;

    tree_id /= base;
  }

  return(edge);
}
