#include <Rcpp/Lightest>
using namespace Rcpp;

#include <algorithm> /* for fill */
#include <string> /* for string (label keys) */
#include <unordered_map> /* for unordered_map */
#include <vector> /* for vector */

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

// Renumber every tree in an UNLABELLED forest so its tips follow `target`'s
// order, entirely in C++. Unlike renumber_tips_batch (one shared permutation
// for a labelled multiPhylo), each tree may carry its own tip.label order, so a
// per-tree permutation is derived from a single hash of `target`. This is the
// fast path for RenumberTips on a plain list / unlabelled multiPhylo, replacing
// a per-tree R loop.
//
// Returns the relabelled list, OR R_NilValue to tell the R caller to fall back
// to the per-tree R path (RenumberTips.phylo). NULL is returned for anything
// this fast path does not handle identically to RenumberTips.phylo: an element
// that is not a phylo-like list or lacks `edge`/`tip.label`; a tree whose label
// set differs from `target` (different length, an unknown or NA label, or a
// non-bijection / duplicate); or duplicate/NA labels in `target` itself. The
// common case — every tree is a permutation of `target` — is handled here.
//
// A tree already in `target` order is returned unchanged (no clone), matching
// RenumberTips.phylo's no-op return.
// [[Rcpp::export]]
SEXP renumber_tips_to(Rcpp::List trees, const Rcpp::CharacterVector target) {
  const int n_tip = target.size();
  const int n_trees = trees.size();
  if (n_tip == 0 || n_trees == 0) return R_NilValue;

  // target label -> 1-based new tip index
  std::unordered_map<std::string, int> pos;
  pos.reserve(n_tip * 2);
  for (int i = 0; i < n_tip; ++i) {
    SEXP s = STRING_ELT(target, i);
    if (s == NA_STRING) return R_NilValue;
    pos.emplace(std::string(CHAR(s)), i + 1);
  }
  if (static_cast<int>(pos.size()) != n_tip) return R_NilValue; // dup target labels

  Rcpp::List result(n_trees);
  std::vector<int> perm(n_tip + 1);
  std::vector<char> used(n_tip + 1);

  for (int t = 0; t < n_trees; ++t) {
    SEXP orig = trees[t];
    if (!Rf_isVectorList(orig)) return R_NilValue;
    Rcpp::List tree_i(orig);
    if (!tree_i.containsElementNamed("tip.label") ||
        !tree_i.containsElementNamed("edge")) return R_NilValue;

    SEXP lab = tree_i["tip.label"];
    if (TYPEOF(lab) != STRSXP || Rf_length(lab) != n_tip) return R_NilValue;

    std::fill(used.begin(), used.end(), char(0));
    bool already = true;
    for (int i = 0; i < n_tip; ++i) {
      SEXP s = STRING_ELT(lab, i);
      if (s == NA_STRING) return R_NilValue;
      auto it = pos.find(std::string(CHAR(s)));
      if (it == pos.end()) return R_NilValue;          // label not in target
      const int new_idx = it->second;
      if (used[new_idx]) return R_NilValue;            // non-bijection (dup label)
      used[new_idx] = 1;
      perm[i + 1] = new_idx;
      if (new_idx != i + 1) already = false;
    }

    if (already) {                                     // no relabel needed
      result[t] = orig;
      continue;
    }

    SEXP orig_class = Rf_getAttrib(orig, R_ClassSymbol);
    Rcpp::List tree_c = Rcpp::clone(tree_i);
    Rcpp::IntegerMatrix edge =
      Rcpp::clone(Rcpp::as<Rcpp::IntegerMatrix>(tree_c["edge"]));
    const int n_edge = edge.nrow();
    for (int j = 0; j < n_edge; ++j) {
      int& child = edge(j, 1);
      if (child <= n_tip) child = perm[child];
    }
    tree_c["edge"] = edge;
    tree_c["tip.label"] = target;
    if (tree_c.hasAttribute("order")) {
      Rcpp::CharacterVector ord = tree_c.attr("order");
      if (ord.size() && ord[0] == "preorder") {
        tree_c.attr("order") = Rcpp::CharacterVector::create("cladewise");
      }
    }
    if (orig_class != R_NilValue) {
      Rf_setAttrib(Rcpp::wrap(tree_c), R_ClassSymbol, orig_class);
    }
    result[t] = tree_c;
  }
  // Preserve list names (the R lapply path this replaces keeps them).
  result.attr("names") = trees.attr("names");
  return result;
}
