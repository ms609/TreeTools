#include <Rcpp.h>
#include <stdexcept>
using namespace Rcpp;
// #define TTDEBUG

#define RETAIN 999

// entry 0 of keep is TRUE if leaf "1" should be retained, false otherwise.
// [[Rcpp::export]]
IntegerMatrix keep_tip (const IntegerMatrix edge, const LogicalVector keep) {
  
  if (edge.ncol() != 2) {
    throw std::invalid_argument("edge must have two columns");
  }
  
  const int
    start_edge = edge.nrow(),
    n_tip = keep.length(),
    all_nodes = start_edge + 2
  ;

  auto
    n_child = std::make_unique<int[]>(all_nodes),
    one_child = std::make_unique<int[]>(all_nodes),
    new_no = std::make_unique<int[]>(all_nodes)
  ;
  
  int next_no = 0;
  for (int i = 0; i != n_tip; ++i) {
    if (keep[i]) {
      n_child[i + 1] = RETAIN;
      new_no[i + 1] = ++next_no;
#ifdef TTDEBUG
      Rcout << " [" << (i + 1) << " -> " << new_no[i + 1] <<"]\n";
#endif
    }
  }
  
  IntegerVector new_child = edge(_, 1);
  int kept_edges = 0;
  for (int i = start_edge; i--; ) {
    const int
      edge_child = edge(i, 1),
      edge_children = n_child[edge_child]
    ;
    if (edge_children) {
      const int edge_parent = edge(i, 0);
      ++n_child[edge_parent];
      if (edge_children == 1) {
        new_child[i] = one_child[edge_child];
        one_child[edge_parent] = one_child[edge_child];
      } else {
        one_child[edge_parent] = edge_child;
        ++kept_edges;
      }
    }
  }
  
#ifdef TTDEBUG
  for (int i = 0; i != start_edge; ++i) {
    Rcout << " - " << (1+i) << ": " << n_child[edge(i, 1)] <<".\n";
  }
#endif
  
  int writing_edge = -1;
  IntegerMatrix ret(kept_edges, 2);
  for (int i = 0; i != start_edge; ++i) {
    const int
      parent = edge(i, 0),
      n_children = n_child[parent]
    ;
    if (n_children) {
      const int child = new_child[i];
      if (n_children > 1) {
        // Record this edge:
        ++writing_edge;
        if (!new_no[parent]) {
          new_no[parent] = ++next_no;
        }
        ret(writing_edge, 0) = new_no[parent];
        if (!new_no[child]) {
          new_no[child] = ++next_no;
        }
        ret(writing_edge, 1) = new_no[child];
#ifdef TTDEBUG
        Rcout << " > Translate: " << parent << "-" << child
              << " --> " << new_no[parent] << "-" << new_no[child] << "\n";
#endif
      }
    }
  }
  
  return ret;
}
