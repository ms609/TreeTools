#include <Rcpp/Lightest>
#include <stdexcept> /* for errors */
#include "../inst/include/TreeTools.h"
using namespace Rcpp;

const uint_fast32_t powers_of_two[32] = {
  0x1, 0x2, 0x4, 0x8,
  0x10, 0x20, 0x40, 0x80,
  0x100, 0x200, 0x400, 0x800,
  0x1000, 0x2000, 0x4000, 0x8000,
  0x10000, 0x20000, 0x40000, 0x80000,
  0x100000, 0x200000, 0x400000, 0x800000,
  0x1000000, 0x2000000, 0x4000000, 0x8000000,
  0x10000000, 0x20000000, 0x40000000, 0x80000000
};

uint_fast32_t bitcounts[65536]; // the bytes representing bit count of each number 0-65535
__attribute__((constructor))
  void initialize_bitcounts() {
    for (int_fast32_t i = 0; i < 65536; i++) {
      int_fast32_t n_bits = 0;
      for (int_fast8_t j = 0; j != 16; j++) {
        if ((i & powers_of_two[j])) ++n_bits;
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
  if (n_tip < 1) {
    Rcpp::stop("nTip < 1");
  }
  if (n_bin != splits.ncol()) {
    Rcpp::stop("nTip does not match split size");
  }

  IntegerVector ret(n_split);
  for (int32 i = n_split; i--; ) {
    for (int32 bin = n_bin; bin--; ) {
      ret[i] += bitcounts[splits(i, bin)];
    }
  }

  return ret;
}
