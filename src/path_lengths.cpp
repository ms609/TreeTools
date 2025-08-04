#include <Rcpp/Lightest>
#include <memory> /* for make_unique */
#include "../inst/include/TreeTools/types.h"
using namespace Rcpp;

#define PARENT(i) edge(i, 0)
#define CHILD(i) edge(i, 1)

#define RTOC(i) (i - 1)
#define CTOR(i) (i + 1)

// edge must be a two-column edge matrix in preorder
// (although seemingly the preorder requirement is not strict;
// see expect_equal(PathLengths(Postorder(bal9)), PathLengths(bal9)))
// [[Rcpp::export]]
NumericMatrix path_lengths(const IntegerMatrix edge, const DoubleVector weight,
                           const LogicalVector init_nas) {
  
  const intx root_node = edge[0];
  const intx n_tip = root_node - 1;
  const intx n_edge = edge.nrow();
  const intx n_vert = n_edge + 1;
  const intx data_dim = n_vert;
  
  Rcpp::NumericMatrix ret = Rcpp::NumericMatrix(Rcpp::no_init(n_vert, n_vert));
  if (init_nas[0]) {
    ret.fill(Rcpp::NumericVector::get_na());
  }
  
  auto parent_of = std::make_unique<intx[]>(n_vert);
  for (intx i = 0; i < n_edge; ++i) {
    const int child_i = RTOC(CHILD(i));
    const int parent_i = PARENT(i);
    parent_of[child_i] = parent_i;
    ret[child_i * data_dim + RTOC(parent_i)] = weight[i];
  }
  
  auto this_path = std::make_unique<intx[]>(n_tip);
  for (intx tip = 1; tip <= n_tip; ++tip) {
    this_path[0] = RTOC(tip);
    intx path_len = 1;
    for(;;) {
      intx this_parent = parent_of[this_path[path_len - 1]];
      if (this_parent) {
        this_path[path_len] = RTOC(this_parent);
      } else {
        break;
      }
      ++path_len;
    }
    
    // span = number of nodes spanned; i.e. edges included - 1
    for (intx span = 1; span < path_len - 1; ++span) {
      for (intx i = 0; i < path_len - span - 1; ++i) {
        const intx* path_i = this_path.get() + i;
        
        const intx start = path_i[span + 1];
        const intx add_to = path_i[span];
        const intx end = path_i[0];
        
        double* ret_ptr = ret.begin();
        
        double* row_end_ptr    = ret_ptr + end    * data_dim;
        double* row_add_to_ptr = ret_ptr + add_to * data_dim;
        
        row_end_ptr[start] = row_add_to_ptr[start] + row_end_ptr[add_to];
      }
    }
  }
  
  return ret;
}
