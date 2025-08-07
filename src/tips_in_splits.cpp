#include <Rcpp/Lightest>

#include "../inst/include/TreeTools.h"
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector tips_in_splits(RawMatrix splits) {
  
  const int32 n_tip = splits.attr("nTip");
  const int32 n_split = splits.nrow();
  const int32 n_bin = (n_tip % 8 == 0 ? 0 : 1) + (n_tip / 8);
  
  if (n_tip < 0) {
    Rcpp::stop("nTip < 0");
  }
  
  if (n_bin != splits.ncol()) {
    Rcpp::stop("nTip does not match split size");
  }

  IntegerVector ret(n_split);
  for (int32 i = 0; i < n_split; ++i) {
    for (int32 bin = 0; bin < n_bin; ++bin) {
      ret[i] += static_cast<int>(__builtin_popcount(splits(i, bin)));
    }
  }

  return ret;
}
