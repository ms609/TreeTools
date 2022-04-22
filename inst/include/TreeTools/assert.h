#ifndef TreeTools_assert_
#define TreeTools_assert_

#ifdef DEBUG
#define ASSERT(x) if (!(x)) {                                  \
Rcpp::Rcerr << "Failed assertion: ";                           \
Rcpp::stop(#x);                                                \
}
#else
#define ASSERT(x) ((void)0)
#endif

#endif
