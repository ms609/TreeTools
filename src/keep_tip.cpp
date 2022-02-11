#include <Rcpp/Lightest>
#include <stdexcept>
using namespace Rcpp;
// #define TTDEBUG

#define RETAIN 9000

#define GET_NEW_NO(n) if (!new_no[n]) new_no[n] = ++next_no;

#define SKIP_EDGE GET_NEW_NO(parent); new_no[child] = new_no[parent]

// entry 0 of keep is TRUE if leaf "1" should be retained, false otherwise.
// [[Rcpp::export]]
IntegerMatrix keep_tip (const IntegerMatrix edge, const LogicalVector keep) {
  if (edge.ncol() != 2) {
    Rcpp::stop("`edge` must have two columns");
  } else {
    IntegerMatrix ret(0, 2);
    if (keep[0]) {
      ret(0, 1) = edge(0, 1);
    }
    return ret;
  }
}
