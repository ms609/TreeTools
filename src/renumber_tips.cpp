#include <Rcpp/Lightest>
using namespace Rcpp;

// Apply a precomputed tip permutation to every tree in batch.
// Returns a list of modified phylo objects (cloned, with new
// edge matrices, updated tip.label, and "preorder" downgraded to "cladewise").
// Preserves the class attribute (typically "phylo") on each element.
// [[Rcpp::export]]
Rcpp::List renumber_tips_batch(
    Rcpp::List trees,
    const Rcpp::IntegerVector perm,
    int n_tip,
    const Rcpp::CharacterVector new_labels
) {
  const int n_trees = trees.size();
  Rcpp::List result(n_trees);

  for (int i = 0; i < n_trees; ++i) {
    SEXP orig = trees[i];
    SEXP orig_class = Rf_getAttrib(orig, R_ClassSymbol);

    // Clone the phylo list (as<List> may strip class; restored below)
    Rcpp::List tree_i = Rcpp::clone(Rcpp::as<Rcpp::List>(orig));

    // Clone and permute the edge matrix
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
    tree_i["edge"] = edge;

    // Replace tip labels (shared across all output trees)
    tree_i["tip.label"] = new_labels;

    // Downgrade "preorder" to "cladewise"
    if (tree_i.hasAttribute("order")) {
      Rcpp::CharacterVector ord = tree_i.attr("order");
      if (ord[0] == "preorder") {
        tree_i.attr("order") = Rcpp::CharacterVector::create("cladewise");
      }
    }

    // Restore class last, after all modifications, to ensure it sticks.
    // Use Rf_setAttrib on the underlying SEXP to bypass Rcpp wrappers.
    if (orig_class != R_NilValue) {
      Rf_setAttrib(Rcpp::wrap(tree_i), R_ClassSymbol, orig_class);
    }

    result[i] = tree_i;
  }

  return result;
}
