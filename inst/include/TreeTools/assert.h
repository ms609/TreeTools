#ifndef TreeTools_assert_
#define TreeTools_assert_

#ifdef DEBUG
#define ASSERT(x) if (!(x)) {                                  \
Rcpp::stop("Assertion failed: " #x);                           \
}
#else
#define ASSERT(x) ((void)0)
#endif

#endif
