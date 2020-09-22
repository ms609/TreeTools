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
  // edge must be BINARY
  // edge must be in preorder
  // #TODO establish the extent to which this really outperforms root_on_node.
  // # Plan to replace with that to reduce future maintenance burden.
  //  [[Rcpp::export]]
  inline IntegerMatrix root_binary(const IntegerMatrix edge, const int outgroup) {

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

  // #TODO Write test cases
  // NB: If specifying internal node by number, note that node numbers will
  // change if tree is not already in preorder.
  //
  //  [[Rcpp::export]]
  inline List root_on_node(const List phy, const int outgroup) {

    IntegerMatrix edge = phy["edge"];

    const intx n_edge = edge.nrow(),
      n_node = phy["Nnode"],
      max_node = n_edge + 1,
      n_tip = max_node - n_node,
      root_node = n_tip + 1
    ;

    edge = preorder_edges_and_nodes(edge(_, 0), edge(_, 1));
    List ret = clone(phy);
    ret.attr("order") = "preorder";
    if (outgroup < 1) throw std::range_error("`outgroup` must be a positive integer");
    if (outgroup > max_node) throw std::range_error("`outgroup` exceeds number of nodes");
    if (outgroup == root_node) {
      ret["edge"] = edge;
      return ret;
    }


    intx* edge_above = new intx[max_node + 1];
    intx root_edges[] = {0, 0};
    intx root_edges_found = 0;

    for (intx i = n_edge; i--; ) {
      edge_above[edge(i, 1)] = i;
      if (edge(i, 0) == root_node) {
        if (root_edges_found < 2) root_edges[root_edges_found] = i;
        ++root_edges_found;
      }
    }

    intx invert_next = edge_above[outgroup];

    if (root_edges_found == 2) { // Root node is vapour, and can be repurposed

      if (edge(root_edges[0], 1) == outgroup ||
          edge(root_edges[1], 1) == outgroup) return phy;
      // #TODO work in situ without clone
      IntegerMatrix new_edge = clone(edge);

      // We'll later add an edge from the now-unallocated root node to the outgroup.
      new_edge(invert_next, 0) = root_node;
      new_edge(invert_next, 1) = edge(invert_next, 0);

      do {
        invert_next = edge_above[edge(invert_next, 0)];
        new_edge(invert_next, 0) = edge(invert_next, 1);
        new_edge(invert_next, 1) = edge(invert_next, 0);
      } while (edge(invert_next, 0) != root_node);

      delete[] edge_above;

      // further root edges must be replaced with root -> outgroup.
      intx spare_edge = (new_edge(root_edges[0], 0) == root_node ? 0 : 1);
      new_edge(invert_next, 1) = edge(root_edges[spare_edge], 1);
      new_edge(root_edges[spare_edge], 1) = outgroup;
      ret["edge"] = preorder_edges_and_nodes(new_edge(_, 0), new_edge(_, 1));

    } else { // Root node will be retained; we need a new root edge

      IntegerMatrix new_edge(n_edge + 1, 2);
      for (int i = n_edge; i--; ) {
        new_edge(i, 0) = edge(i, 0);
        new_edge(i, 1) = edge(i, 1);
      }
      const intx new_root = max_node + 1;
      new_edge(n_edge, 0) = new_root;
      new_edge(n_edge, 1) = outgroup;

      new_edge(invert_next, 0) = new_root;
      new_edge(invert_next, 1) = edge(invert_next, 0);

      while (edge(invert_next, 0) != root_node) {
        invert_next = edge_above[edge(invert_next, 0)];
        new_edge(invert_next, 0) = edge(invert_next, 1);
        new_edge(invert_next, 1) = edge(invert_next, 0);
      }
      delete[] edge_above;

      ret["Nnode"] = n_node + 1;
      ret["edge"] = preorder_edges_and_nodes(new_edge(_, 0), new_edge(_, 1));

    }
    // #TODO there is probably a clever way to avoid doing a full preorder rewriting.
    return ret;
  }
}

#endif
