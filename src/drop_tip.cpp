#include <Rcpp.h>
#include <stdexcept>
using namespace Rcpp;

// Docs do not assert that edge must be in preorder, but warn against
// "unconventional" ordering.  (See `nasty` tree in test cases.)
// We require that the root node is numbered nTip + 1.
// drop is a vector of integers between 1 and nTip marking tips to drop.
// [[Rcpp::export]]
IntegerMatrix drop_tip (const IntegerMatrix edge, const IntegerVector drop) {

  if (edge.ncol() != 2) {
    throw std::invalid_argument("edge must have two columns");
  }
  const int
    start_edge = edge.nrow(),
    all_nodes = start_edge + 2
  ;

  LogicalVector
    dropped_edge(start_edge),
    dropped_node(all_nodes)
  ;

  IntegerVector new_child = edge(_, 1);
  auto
    parent_of = std::make_unique<int[]>(all_nodes),
    n_children = std::make_unique<int[]>(all_nodes),
    edge_above = std::make_unique<int[]>(all_nodes)
  ;
  for (int i = start_edge; i--; ) {
    const int
      parent = edge(i, 0),
      child = new_child[i]
    ;
    parent_of[child] = parent;
    ++n_children[parent];
    edge_above[child] = i;
  }
  int root = start_edge / 2 + 2; // Minimum possible value
  while (edge_above[root]) {
    ++root;
  }
  const bool rooted = n_children[root] == 2;

  for (int i = drop.length(); i--; ) {
    const int tip = drop(i);
    dropped_node[tip] = true;
    dropped_edge[edge_above[tip]] = true;

    const int node_parent = parent_of[tip];
#ifdef MSDEBUG
    Rcout << " Dropping tip " << tip << " on edge 1+" << edge_above[tip] << ".\n";
    Rcout << "  Parent of tip = " << node_parent << ", instance "
          << n_children[node_parent] <<"\n";
#endif
    parent_of[tip] = 0;
    if (--n_children[node_parent] == 1) {
      int child = 0;
      do {++child;} while (parent_of[child] != node_parent);
      dropped_edge[edge_above[child]] = true;
      const int
        lost_node = parent_of[child],
        child_destination = edge_above[lost_node]
      ;
#ifdef MSDEBUG
      Rcout << " Lost node " << lost_node << "; moving child " << child
            << " to edge 1+"
            << edge_above[lost_node] << " and setting parent to "
            << parent_of[lost_node] << ".\n";
#endif
      dropped_node[lost_node] = true;
      edge_above[child] = child_destination;
      new_child[child_destination] = child;
      parent_of[child] = parent_of[lost_node];
      parent_of[lost_node] = 0;
    }
  }


  int dropped = 0;
  auto new_no = std::make_unique<int[]>(all_nodes);
  for (int i = 1; i != all_nodes; i++) {
    if (dropped_node[i]) {
      ++dropped;

#ifdef MSDEBUG
      Rcout << " Summary: Dropped another node: " << i <<". That's "
            << dropped << " altogether.\n";
#endif
#ifndef NDEBUG
      new_no[i] = -i;
#endif

    } else {
      new_no[i] = i - dropped;

#ifdef MSDEBUG
      Rcout << "  New no for " << i <<": " << new_no[i] << " .\n";
#endif
    }
  }

  const int new_root = root - drop.length();
  if (!rooted && new_root > 3) {
    // Check to see whether we've inadvertently been left with a root node
    // Example: drop tip 1 from "(1, 2, (3, 4));" -> "(2, (3, 4));"
    int
      root_order = 0,
      collapse
    ;
    for (int i = start_edge; i--; ) {
      if (new_no[edge(i, 0)] == new_root) {
        ++root_order;
        if (new_no[edge(i, 1)] > new_root) {
          collapse = i;
        }
      }
    }
    if (root_order == 2) {
#ifdef MSDEBUG
      Rcout << "  Removing root node by removing edge " << collapse 
            << ": " << new_no[edge(collapse, 0)] 
            << " - " << new_no[edge(collapse, 1)] << ".\n";
#endif
      dropped_edge[collapse] = true;
      ++dropped;
      new_no[edge(collapse, 1)] = new_root;
    }
  }
  
  IntegerMatrix ret(start_edge - dropped, 2);
  for (int i = start_edge; i--; ) {
    if (dropped_edge[i]) {
      --dropped;
    } else {
      assert(new_no[edge(i, 0)] > 0);
      assert(new_no[new_child[i]] > 0);
      ret(i - dropped, 0) = new_no[edge(i, 0)];
      ret(i - dropped, 1) = new_no[new_child[i]];
    }
  }

  return ret;

}
