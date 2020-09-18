#ifndef TreeTools_root_tree_
#define TreeTools_root_tree_

#include <Rcpp.h>
#include "renumber_tree.h"
#include "types.h"

namespace TreeTools {
  using namespace Rcpp;
  extern inline IntegerMatrix preorder_edges_and_nodes(const IntegerVector parent,
                                                       const IntegerVector child);

  // #TODO Write test cases
  // edge must be in preorder
  //  [[Rcpp::export]]
  inline IntegerMatrix root_on_node(const IntegerMatrix edge, const int outgroup) {

    if (edge(0, 1) == outgroup) return edge;

    const intx n_edge = edge.nrow(),
      n_node = n_edge / 2,
      n_tip = n_node + 1,
      root_node = n_tip + 1,
      max_node = n_node + n_tip;

    if (outgroup < 1) throw std::range_error("`outgroup` must be a positive integer");
    if (outgroup > max_node) throw std::range_error("`outgroup` exceeds number of nodes");
    if (outgroup == root_node) return edge;

    intx* edge_above = new intx[max_node + 1];
    intx root_edges[2] = {0, 0};

    for (intx i = n_edge; i--; ) {

      edge_above[edge(i, 1)] = i;

      if (edge(i, 0) == root_node) {
        if (edge(i, 1) == outgroup) {
          delete[] edge_above;
          return edge;
        }
        root_edges[root_edges[1] ? 0 : 1] = i;
      }

    }

    IntegerMatrix ret = clone(edge);
    intx invert_next = edge_above[outgroup];

    // We'll later add an edge from the now-unallocated root node to the outgroup.
    ret(invert_next, 0) = root_node;
    ret(invert_next, 1) = edge(invert_next, 0);

    do {
      invert_next = edge_above[edge(invert_next, 0)];
      ret(invert_next, 0) = edge(invert_next, 1);
      ret(invert_next, 1) = edge(invert_next, 0);
    } while (edge(invert_next, 0) != root_node);

    delete[] edge_above;

    // second root i.e. 16 -- 24 must be replaced with root -> outgroup.
    intx spare_edge = (ret(root_edges[0], 0) == root_node ? 0 : 1);
    ret(invert_next, 1) = edge(root_edges[spare_edge], 1);
    ret(root_edges[spare_edge], 1) = outgroup;

    return preorder_edges_and_nodes(ret(_, 0), ret(_, 1));
  }
}

#endif
