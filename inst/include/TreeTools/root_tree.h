#ifndef TreeTools_root_tree_
#define TreeTools_root_tree_

#include <Rcpp/Lightest>
#include <memory> /* for std::unique_ptr, make_unique */
#include <stdexcept> /* for errors */
#include "assert.h" /* for ASSERT */
#include "renumber_tree.h"
#include "types.h"

namespace TreeTools {
  extern inline Rcpp::IntegerMatrix preorder_edges_and_nodes(
      const Rcpp::IntegerVector parent,
      const Rcpp::IntegerVector child);

  extern inline Rcpp::List preorder_weighted(
    const Rcpp::IntegerVector parent,
    const Rcpp::IntegerVector child,
    const Rcpp::DoubleVector weight);

  // #TODO Write test cases
  // edge must be BINARY
  // edge must be in preorder
  // #TODO establish the extent to which this really outperforms root_on_node.
  // # Plan to replace with that to reduce future maintenance burden.
  //  [[Rcpp::export]]
  inline Rcpp::IntegerMatrix root_binary(const Rcpp::IntegerMatrix edge,
                                         const int outgroup) {

    if (edge(0, 1) == outgroup) return edge;

    const intx n_edge = edge.nrow(),
      n_node = n_edge / 2,
      n_tip = n_node + 1,
      root_node = n_tip + 1,
      max_node = n_node + n_tip;

    if (outgroup < 1) {
      Rcpp::stop("`outgroup` must be a positive integer");
    }
    if (outgroup > max_node) {
      Rcpp::stop("`outgroup` exceeds number of nodes");
    }
    if (outgroup == root_node) {
      return edge;
    }
    

    std::unique_ptr<intx[]> edge_above = std::make_unique<intx[]>(max_node + 1);
    intx root_edges[2] = {0, 0};

    for (intx i = n_edge; i--; ) {

      edge_above[edge(i, 1)] = i;

      if (edge(i, 0) == root_node) {
        if (edge(i, 1) == outgroup) {
          return edge;
        }
        root_edges[root_edges[1] ? 0 : 1] = i;
      }

    }

    Rcpp::IntegerMatrix ret = Rcpp::clone(edge);
    intx invert_next = edge_above[outgroup];

    // We'll later add an edge from the now-unallocated root node to the outgroup.
    ret(invert_next, 0) = root_node;
    ret(invert_next, 1) = edge(invert_next, 0);

    do {
      invert_next = edge_above[edge(invert_next, 0)];
      ret(invert_next, 0) = edge(invert_next, 1);
      ret(invert_next, 1) = edge(invert_next, 0);
    } while (edge(invert_next, 0) != root_node);

    // second root i.e. 16 -- 24 must be replaced with root -> outgroup.
    intx spare_edge = (ret(root_edges[0], 0) == root_node ? 0 : 1);
    ret(invert_next, 1) = edge(root_edges[spare_edge], 1);
    ret(root_edges[spare_edge], 1) = outgroup;

    return preorder_edges_and_nodes(ret(Rcpp::_, 0), ret(Rcpp::_, 1));
  }

  // #TODO Write test cases
  // NB: If specifying internal node by number, note that node numbers will
  // change if tree is not already in preorder.
  // NB: root_node must == n_tip + 1
  //
  //  [[Rcpp::export]]
  inline Rcpp::List root_on_node(const Rcpp::List phy, const int outgroup) {

    Rcpp::IntegerMatrix edge = phy["edge"];
    Rcpp::NumericVector weight;

    const intx
      n_edge = edge.nrow(),
      n_node = phy["Nnode"],
      max_node = n_edge + 1,
      n_tip = max_node - n_node,
      root_node = n_tip + 1
    ;
    const bool weighted = phy.containsElementNamed("edge.length");

    if (weighted) {
      Rcpp::List reweighted = preorder_weighted(
        edge(Rcpp::_, 0),
        edge(Rcpp::_, 1),
        phy["edge.length"]
      );
      Rcpp::IntegerMatrix edge = reweighted[0];
      Rcpp::NumericVector weight = reweighted[1];
    } else {
      edge = preorder_edges_and_nodes(edge(Rcpp::_, 0), edge(Rcpp::_, 1));
    }
    if (outgroup < 1) {
      Rcpp::stop("`outgroup` must be a positive integer");
    }
    if (outgroup > max_node) {
      Rcpp::stop("`outgroup` exceeds number of nodes");
    }
    Rcpp::List ret = Rcpp::clone(phy);
    ret.attr("order") = "preorder";
    if (outgroup == root_node) {
      ret["edge"] = edge;
      if (weighted) {
        ret["edge.weight"] = weight;
      }
      return ret;
    }


    auto edge_above = std::make_unique<intx[]>(max_node + 1);
    intx root_edges[] = {0, 0};
    intx root_edges_found = 0;

    for (intx i = n_edge; i--; ) {
      edge_above[edge(i, 1)] = i;
      if (edge(i, 0) == root_node) {
        if (root_edges_found < 2) {
          root_edges[root_edges_found] = i;
        }
        ++root_edges_found;
      }
    }

    intx invert_next = edge_above[outgroup];

    if (root_edges_found == 2) { // Root node is vapour, and can be repurposed

      if (edge(root_edges[0], 1) == outgroup ||
          edge(root_edges[1], 1) == outgroup) {
        return phy;
      }
      // #TODO work in situ without clone?
      Rcpp::IntegerMatrix new_edge = clone(edge);

      // We'll later add an edge from the now-unallocated root node to the outgroup.
      new_edge(invert_next, 0) = root_node;
      new_edge(invert_next, 1) = edge(invert_next, 0);

      do {
        invert_next = edge_above[edge(invert_next, 0)];
        new_edge(invert_next, 0) = edge(invert_next, 1);
        new_edge(invert_next, 1) = edge(invert_next, 0);
      } while (edge(invert_next, 0) != root_node);

      // Further root edges must be replaced with root -> outgroup.
      intx spare_edge = (new_edge(root_edges[0], 0) == root_node ? 0 : 1);
      new_edge(invert_next, 1) = edge(root_edges[spare_edge], 1);
      new_edge(root_edges[spare_edge], 1) = outgroup;
      if (weighted) {
        Rcpp::List preorder_res;
        preorder_res = preorder_weighted(new_edge(Rcpp::_, 0),
                                         new_edge(Rcpp::_, 1),
                                         phy["edge.length"]);
        ret["edge"] = preorder_res[0];
        ret["edge.length"] = preorder_res[1];
      } else {
        ret["edge"] = preorder_edges_and_nodes(new_edge(Rcpp::_, 0),
                                               new_edge(Rcpp::_, 1));
      }

    } else { // Root node will be retained; we need a new root edge

      Rcpp::IntegerMatrix new_edge(n_edge + 1, 2);
      Rcpp::NumericVector new_wt(n_edge + 1);
      if (weighted) {
        weight = phy["edge.length"];
        for (int i = n_edge; i--; ) {
          new_wt[i] = weight[i];
        }
        ASSERT(new_wt(n_edge) == 0);
      }
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

      ret["Nnode"] = n_node + 1;
      if (weighted) {
        Rcpp::List preorder_res;
        preorder_res = preorder_weighted(
          new_edge(Rcpp::_, 0),
          new_edge(Rcpp::_, 1),
          new_wt);
        ret["edge"] = preorder_res[0];
        ret["edge.length"] = preorder_res[1];
      } else {
        ret["edge"] = preorder_edges_and_nodes(new_edge(Rcpp::_, 0),
                                               new_edge(Rcpp::_, 1));
      }
      
    }
    // #TODO there is probably a clever way to avoid doing a full preorder rewriting.
    return ret;
  }
}

#endif
