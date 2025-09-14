#include <Rcpp/Lightest>
#include <vector> /* for errors */
#include <stdexcept> /* for errors */
#include "../inst/include/TreeTools/n_cherries.h"

// [[Rcpp::export]]
Rcpp::IntegerVector n_cherries_wrapper(const Rcpp::IntegerVector parent,
                                       const Rcpp::IntegerVector child,
                                       const int nTip) {
  try {
    const size_t n_edge = parent.size();
    if (child.size() != (int)n_edge) {
      Rcpp::stop("`parent` and `child` must be the same length");
    }
    
    // Call your C++ function
    int result = TreeTools::n_cherries(parent.begin(), child.begin(),
                                       n_edge, nTip);
    
    // Return the result as an R integer
    return Rcpp::wrap(result);
  } catch (const std::exception& e) {
    // Catch any standard exception and throw it as an R error
    Rcpp::stop(e.what());
  } catch (...) {
    // Catch any other exceptions
    Rcpp::stop("An unknown error occurred in n_cherries_wrapper.");
  }
}