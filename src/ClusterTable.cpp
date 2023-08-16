#include "../inst/include/TreeTools/ClusterTable.h"
#include <Rcpp/Lightest>

// Modelled on https://CRAN.R-project.org/package=Rcpp/vignettes/Rcpp-modules.pdf
// [[Rcpp::export]]
SEXP ClusterTable_new(Rcpp::List phylo) {
  Rcpp::XPtr<TreeTools::ClusterTable> 
    ptr(new TreeTools::ClusterTable(phylo), true);
  
  return ptr;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix ClusterTable_matrix(SEXP xp) {
  Rcpp::XPtr<TreeTools::ClusterTable> ptr(xp);
  return ptr->X_contents();
}

// [[Rcpp::export]]
Rcpp::IntegerVector ClusterTable_decode(SEXP xp) {
  Rcpp::XPtr<TreeTools::ClusterTable> ptr(xp);
  return ptr->X_decode();
}
