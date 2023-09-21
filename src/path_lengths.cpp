#include <Rcpp/Lightest>
#include <memory> /* for make_unique */
#include "../inst/include/TreeTools/types.h"
using namespace Rcpp;

#define PARENT(i) edge(i, 0)
#define CHILD(i) edge(i, 1)

// edge must be a two-column edge matrix in preorder
// [[Rcpp::export]]
NumericMatrix path_lengths(const IntegerMatrix edge, const DoubleVector weight) {
  const intx
    root_node = edge[0],
    n_tip = root_node - 1,
    n_edge = edge.nrow(),
    n_vert = n_edge + 1,
    r_to_c = 1
  ;
  
  NumericMatrix ret(n_vert + 1, n_vert + 1);
  ret.fill(NumericVector::get_na());
  auto
    parent_of = std::make_unique<intx[]>(n_vert + r_to_c),
    parent_edge = std::make_unique<intx[]>(n_vert + r_to_c)
  ;
  for (intx i = n_edge; i--; ) {
    parent_of[CHILD(i)] = PARENT(i);
    parent_edge[CHILD(i)] = i;
    ret(PARENT(i), CHILD(i)) = weight[i];
  }
  auto this_path = std::make_unique<intx[]>(n_tip);
  for (intx tip = 1; tip <= n_tip; ++tip) {
    this_path[0] = tip;
    intx path_len = 1;
    for(;;) {
      intx this_parent = parent_of[this_path[path_len - 1]];
      if (this_parent) {
        this_path[path_len] = this_parent;
      } else {
        break;
      }
      ++path_len;
    }
    // span = number of nodes spanned; i.e. edges included - 1
    for (intx span = 1; span < path_len - 1; ++span) {
      for (intx i = 0; i != path_len - span - 1; ++i) {
        const intx
          start = this_path[i + span + 1],
          add_to = this_path[i + span],
          end = this_path[i]
        ;
        ret(start, end) = ret(start, add_to) + ret(add_to, end);
      }
    }
  }
  return ret(Range(1, n_vert), Range(1, n_vert));
}
