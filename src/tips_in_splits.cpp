#include <Rcpp/Lightest>

#include "../inst/include/TreeTools.h"

// [[Rcpp::export]]
Rcpp::IntegerVector tips_in_splits(Rcpp::RawMatrix splits) {
  
  const int32 n_tip = splits.attr("nTip");
  const int32 n_split = splits.nrow();
  const int32 n_bin = (n_tip % 8 == 0 ? 0 : 1) + (n_tip / 8);
  
  if (n_tip < 0) {
    Rcpp::stop("nTip < 0");
  }
  
  if (n_bin != splits.ncol()) {
    Rcpp::stop("nTip does not match split size");
  }

  Rcpp::IntegerVector ret(n_split);
  for (int32 i = 0; i < n_split; ++i) {
    const unsigned char* splits_i = splits.begin() + i; 
    for (int32 bin = 0; bin < n_bin; ++bin) {
      const unsigned char* in_bin = splits_i + bin * n_split;
      ret[i] += static_cast<int>(__builtin_popcount(*in_bin));
    }
  }

  return ret;
}
