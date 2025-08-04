#include <Rcpp/Lightest>

#include <stdexcept> /* for errors */
#include "../inst/include/TreeTools.h"
using namespace Rcpp;

uint_fast32_t bitcounts[65536]; // the bytes representing bit count of each number 0-65535
__attribute__((constructor))
  void initialize_bitcounts() {
    for (int_fast32_t i = 0; i < 65536; i++) {
      int_fast32_t n_bits = 0;
      for (int_fast8_t j = 0; j != 16; j++) {
        if (i & uint_fast32_t(1) << j) ++n_bits;
      }
      bitcounts[i] = n_bits;
    }
  }

// [[Rcpp::export]]
IntegerVector tips_in_splits(RawMatrix splits) {
  const int32
    n_tip = splits.attr("nTip"),
    n_split = splits.nrow(),
    n_bin = (n_tip % 8 == 0 ? 0 : 1) + (n_tip / 8)
  ;
  if (n_tip < 0) {
    Rcpp::stop("nTip < 0");
  }
  if (n_bin != splits.ncol()) {
    Rcpp::stop("nTip does not match split size");
  }

  IntegerVector ret(n_split);
  for (int32 i = n_split; i--; ) {
    for (int32 bin = n_bin; bin--; ) {
      ret[i] += decltype(ret[0])(bitcounts[splits(i, bin)]);
    }
  }

  return ret;
}
