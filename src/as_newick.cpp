#include <Rcpp/Lighter>
#include <string>
#include "../inst/include/TreeTools/renumber_tree.h"

using namespace Rcpp;

const unsigned int MAX_NEWICK_NODES = 8192;

void close_node (const unsigned int top_node,
                 const unsigned int bottom_node,
                 bool is_open[],
                 const unsigned int n_tip,
                 std::string& ret) {
  for (unsigned int j = top_node; j > bottom_node; j--) {
    if (is_open[j - n_tip]) {
      ret.append(")");
      is_open[j - n_tip] = false;
    }
  }
}

// `edge` should be the $edge element of a `phylo` object, minus one.
// [[Rcpp::export]]
CharacterVector as_newick(const IntegerMatrix edge) {
  const unsigned int n_row = edge.nrow();
  if (n_row > MAX_NEWICK_NODES * 2 - 1) {
    Rcpp::stop("Too many nodes for as_newick");
  }
  if (!n_row) {
    return CharacterVector::create(";");
  }
  if (edge.ncol() != 2) {
    Rcpp::stop("`edge` must have two columns");
  }
  if (Rcpp::min(edge) != 0) {
    if (Rcpp::min(edge) == NA_INTEGER) {
      Rcpp::stop("`edge` may not contain NA entries");
    } else {
      Rcpp::stop("`min(edge)` must be zero");
    }
  }
  if (Rcpp::max(edge) != edge.nrow()) {
    Rcpp::stop("`edge` is malformed");
  }

  unsigned int last_node = 0;
  bool is_open[MAX_NEWICK_NODES];
  std::string ret;
  ret.reserve(n_row * 6); /* Four characters for tip + "," + "(" or ")" */
  IntegerMatrix preorder =
    TreeTools::preorder_edges_and_nodes(edge(_, 0) + 1, edge(_, 1) + 1) - 1;
  const unsigned int n_tip = preorder(0, 0);

  for (unsigned int i = 0; i != n_row; i++) {
    unsigned int this_node = preorder(i, 0), this_child = preorder(i, 1);
    if (this_node == last_node) {
      ret.append(",");
    } else if (this_node > last_node) {
      ret.append("(");
      is_open[this_node - n_tip] = true;
    } else {
      close_node(last_node, this_node, is_open, n_tip, ret);
      ret.append(",");
    }
    if (this_child < n_tip) {
      ret.append(std::to_string(this_child));
    }
    last_node = this_node;
  }
  close_node(preorder(preorder.nrow() - 1, 0), n_tip - 1, is_open, n_tip, ret);
  ret.append(";");
  return CharacterVector::create(ret);
}
