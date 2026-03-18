#include <Rcpp/Lightest>
using namespace Rcpp;

// Apply a precomputed tip permutation to every tree's edge matrix in batch.
// perm: 1-indexed permutation vector (matchOrder from R's match())
// n_tip: number of tips
// Returns a list of new (cloned) edge matrices with remapped tip indices.
// [[Rcpp::export]]
Rcpp::List renumber_tips_batch(
    Rcpp::List trees,
    const Rcpp::IntegerVector perm,
    int n_tip
) {
  const int n_trees = trees.size();
  Rcpp::List result(n_trees);

  for (int i = 0; i < n_trees; ++i) {
    Rcpp::List tree_i = Rcpp::as<Rcpp::List>(trees[i]);
    Rcpp::IntegerMatrix edge = Rcpp::clone(
      Rcpp::as<Rcpp::IntegerMatrix>(tree_i["edge"])
    );
    const int n_edge = edge.nrow();

    for (int j = 0; j < n_edge; ++j) {
      int& child = edge(j, 1);
      if (child <= n_tip) {
        child = perm[child - 1];
      }
    }

    result[i] = edge;
  }

  return result;
}
