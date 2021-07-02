#include <Rcpp.h>
using namespace Rcpp;

#include "../inst/include/TreeTools/types.h"
#include "../inst/include/TreeTools/renumber_tree.h"

#define ASSERT_NODE_COUNT(x) assert((x) - start_tip - 1 >= 0); assert((x) < (start_tip + start_node + 1))
#define NODE_COUNT(x) node_count[(x) - start_tip - 1]

#define RAISE_EDGE(i) assert((i) >= n_deleted);                       \
  if (n_deleted) ret((i) - n_deleted, _) = ret((i), _)

// edges must be in preorder.
// drop is a vector of integers between 1 and nTip marking tips to drop.
// [[Rcpp::export]]
IntegerMatrix drop_tip2 (const IntegerMatrix preorder, const IntegerVector drop) {

  const int32
    root_node = preorder(0, 0),
    start_tip = root_node - 1,
    start_edge = preorder.nrow(),
    start_node = start_edge + 1 - start_tip
  ;

  int32
    still_to_drop = drop.length(),
      n_deleted = 0,
      ret_edges = start_edge - still_to_drop,
      next_no = start_tip - still_to_drop + 1
  ;

  IntegerVector droppers = clone(drop);
  IntegerMatrix ret(ret_edges, 2);

  std::unique_ptr<int32[]>
    node_count = std::make_unique<int32[]>(start_node),
      new_no = std::make_unique<int32[]>(start_tip + start_node + 1)
    ;


  std::sort(droppers.begin(), droppers.end());
  for (int32 i = 1; i != start_tip + 1; i++) {
    if (n_deleted != droppers.length() && i == droppers[n_deleted]) {
      ++n_deleted;
    } else {
      new_no[i] = i - n_deleted;
    }
  }


  for (int32 i = start_edge; i--; ) {
    const int32 child = preorder(i, 1);
    assert(child < (start_tip + start_node + 1));
    if (child <= start_tip && !new_no[child]) {
      still_to_drop--;
    } else {
      ASSERT_NODE_COUNT(preorder(i, 0));
      NODE_COUNT(preorder(i, 0))++;
      ret(i - still_to_drop, _) = preorder(i, _);
    }
  }

  do {
    n_deleted = (NODE_COUNT(ret(0, 0)) < 2) ? 1 : 0;
    for (int32 i = n_deleted; i < ret_edges; i++) { // preorder traversal
      const int32 child = ret(i, 1);
      if (child > start_tip) {
        ASSERT_NODE_COUNT(ret(i, 1));
        if (NODE_COUNT(ret(i, 1)) == 0) {
          ASSERT_NODE_COUNT(ret(i, 0));
          NODE_COUNT(ret(i, 0))--;
          ++n_deleted;
        } else if (NODE_COUNT(ret(i, 1)) == 1) {
          NODE_COUNT(ret(i, 1))--;
          if (n_deleted) {
            ret(i - n_deleted, 0) = ret(i, 0);
            ret(i - n_deleted, 1) = ret(i + 1, 1);
          } else {
            ret(i, 1) = ret(i + 1, 1);
          }
          ++n_deleted;
          ++i;
        } else {
          RAISE_EDGE(i);
        }
      } else {
        RAISE_EDGE(i);
      }
    }
    ret_edges -= n_deleted;
  } while (n_deleted);

  new_no[ret(0, 0)] = next_no++;
  for (int32 i = 0; i != ret_edges; i++) {
    const int32 parent = ret(i, 0), child = ret(i, 1);
    assert(parent < (start_tip + start_node + 1));
    assert(child < (start_tip + start_node + 1));
    assert(child > 0);
    assert(new_no[parent]);
    if (!new_no[child]) {
      new_no[child] = next_no++;
    }

    //
    //     Rcout << "    [" << i << ",]  " << new_no[parent] << "[" << parent << "] "
    //           << new_no[child] << "[" << child << "]\n";

    ret(i, 0) = new_no[parent];
    ret(i, 1) = new_no[child];
  }


  return ret(Range(0, ret_edges - 1), _);

}

// edges must be in preorder.
// drop is a vector of integers between 1 and nTip marking tips to drop.
// [[Rcpp::export]]
IntegerMatrix drop_tip (const IntegerMatrix edge, const IntegerVector drop) {

  const int32
    start_edge = edge.nrow(),
    all_nodes = start_edge + 2
  ;

  std::unique_ptr<int32[]>
    parent_of = std::make_unique<int32[]>(all_nodes),
    as_parent = std::make_unique<int32[]>(all_nodes),
    child_on =  std::make_unique<int32[]>(all_nodes),
    new_child = std::make_unique<int32[]>(start_edge),
    new_no = std::make_unique<int32[]>(all_nodes)
  ;
  std::unique_ptr<bool[]>
    dropped_edge = std::make_unique<bool[]>(start_edge),
    dropped_node = std::make_unique<bool[]>(all_nodes)
  ;


  for (int32 i = start_edge; i--; ) {
    const int32 parent = edge(i, 0), child = edge(i, 1);
    parent_of[child] = parent;
    ++as_parent[parent];
    child_on[child] = i;
    new_child[i] = child;
  }

  for (int32 i = drop.length(); i--; ) {
    const int32 tip = drop(i);
    dropped_node[tip] = true;
    dropped_edge[child_on[tip]] = true;

    const int32 node_parent = parent_of[tip];
#ifdef MSDEBUG
    Rcout << " Dropping tip " << tip << " on edge 1+"<<child_on[tip] <<".\n";
    Rcout << "  Parent of tip = " << node_parent << ", instance "
          << as_parent[node_parent] <<"\n";
#endif
    parent_of[tip] = 0;
    if (--as_parent[node_parent] == 1) {
      int32 child = 0;
      do {++child;} while (parent_of[child] != node_parent);
      dropped_edge[child_on[child]] = true;
      const int32
        lost_node = parent_of[child],
        child_destination = child_on[lost_node]
      ;
#ifdef MSDEBUG
      Rcout << " Lost node " << lost_node << "; moving child " << child
            << " to edge 1+"
            << child_on[lost_node] << " and setting parent to "
            << parent_of[lost_node] << ".\n";
#endif
      dropped_node[lost_node] = true;
      child_on[child] = child_destination;
      new_child[child_destination] = child;
      parent_of[child] = parent_of[lost_node];
      parent_of[lost_node] = 0;
    }
  }


  int32 dropped = 0;
  for (int32 i = 1; i != all_nodes; i++) {
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

  IntegerMatrix ret(start_edge - dropped, 2);
  for (int32 i = start_edge; i--; ) {
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
