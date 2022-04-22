#include <Rcpp/Lightest>
#include <string>
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

// Edge should be the $edge element of a `phylo` object, minus one.
// [[Rcpp::export]]
CharacterVector as_newick(IntegerMatrix edge) {
  if ((unsigned int) edge.nrow() > MAX_NEWICK_NODES * 2 - 1) {
    Rcpp::stop("Too many nodes for as_newick");
  }

  const unsigned int n_tip = edge(0, 0);
  unsigned int last_node = 0;
  bool is_open[MAX_NEWICK_NODES];
  std::string ret;
  ret.reserve(edge.nrow() * 6); /* Four characters for tip + "," + "(" or ")" */

  for (unsigned int i = 0; i < (unsigned int) edge.nrow(); i++) {
    unsigned int this_node = edge(i, 0), this_child = edge(i, 1);
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
  close_node(edge(edge.nrow() - 1, 0), n_tip - 1, is_open, n_tip, ret);
  ret.append(";");
  return CharacterVector::create(ret);
}
