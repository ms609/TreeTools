#include <Rcpp/Lightest>
#include <memory> /* for make_unique */
#include "../inst/include/TreeTools/types.h"
using namespace Rcpp;

#define PARENT(i) edge(i, 0)
#define CHILD(i) edge(i, 1)

// edge must be a two-column edge matrix in preorder
// [[Rcpp::export]]
NumericMatrix path_lengths(const IntegerMatrix edge, const DoubleVector weight) {
  
  const intx root_node = edge[0];
  const intx n_tip = root_node - 1;
  const intx n_edge = edge.nrow();
  const intx n_vert = n_edge + 1;
  const intx data_dim = n_vert + 1;
  constexpr intx r_to_c = 1;
  
  std::vector<double> data(data_dim * data_dim, NumericVector::get_na());
  
  auto parent_of = std::make_unique<intx[]>(n_vert + r_to_c);
  auto parent_edge = std::make_unique<intx[]>(n_vert + r_to_c);
  for (intx i = n_edge; i--; ) {
    parent_of[CHILD(i)] = PARENT(i);
    parent_edge[CHILD(i)] = i;
    data[CHILD(i) * data_dim + PARENT(i)] = weight[i];
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
        const intx start = this_path[i + span + 1];
        const intx add_to = this_path[i + span];
        const intx end = this_path[i];
        
        data[end * data_dim + start] = data[add_to * data_dim + start] + 
          data[end * data_dim + add_to];
      }
    }
  }
  NumericMatrix ret(data_dim, data_dim, data.begin());
  return ret(Range(1, n_vert), Range(1, n_vert));
}
