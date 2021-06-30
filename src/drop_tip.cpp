#include <Rcpp.h>
using namespace Rcpp;

#include "../inst/include/TreeTools/types.h"
#include "../inst/include/TreeTools/renumber_tree.h"

#define ASSERT_NODE_COUNT(x) assert((x) - start_tip - 1 >= 0); assert((x) < (start_tip + start_node + 1))
#define NODE_COUNT(x) node_count[(x) - start_tip - 1]

#define RAISE_EDGE(i) assert((i) >= n_deleted);                       \
  Rcout << " Raising edge " << i << " to " << (i - n_deleted) << ": " \
        << int(ret(i, 0)) << " " << int(ret(i, 1)) << "\n";           \
  if (n_deleted) ret((i) - n_deleted, _) = ret((i), _)

// [[Rcpp::export]]
IntegerMatrix drop_tip (const IntegerMatrix edge, const IntegerVector drop) {
  const IntegerMatrix preorder = TreeTools::preorder_edges_and_nodes(edge(_, 0), edge(_, 1));

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
    next_no = start_tip + 2 - still_to_drop // tips + root node
  ;

  IntegerVector droppers = clone(drop);
  IntegerMatrix ret(ret_edges, 2);

  std::unique_ptr<int32[]>
    node_count = std::make_unique<int32[]>(start_node),
  ;

  Rcout << "Starting on tree with " << start_tip << " tips, root = "
        << root_node <<", " << start_node << " nodes, "
        << start_edge << " edges.\n";

  std::sort(droppers.begin(), droppers.end());
  for (int32 i = 1; i != start_tip + 2; i++) { // tips + root node
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

  for (int i = 0; i != start_node; i++) {
    Rcout << node_count[i] << " instances of node ";
    Rcout << (i + start_tip + 1) << "\n";
  }

  do {
    n_deleted = 0;
    Rcout << "\n" << ret_edges << " edges in ret\n\n";
    for (int32 i = 0; i < ret_edges; i++) { // postorder traversal
      Rcout << "   [" << i << ",]  " << ret(i, 0) << "   " << ret(i, 1) << "\n";
      const int32 child = ret(i, 1);
      if (child > start_tip) {
        ASSERT_NODE_COUNT(ret(i, 1));
        if (NODE_COUNT(ret(i, 1)) == 0) {
          Rcout << "Deleting dead end edge " << i <<"\n";
          ASSERT_NODE_COUNT(ret(i, 0));
          NODE_COUNT(ret(i, 0))--;
          ++n_deleted;
        } else if (NODE_COUNT(ret(i, 1)) == 1) {
          Rcout << "Deleting singleton edge " << i << ": "
          << ret(i, 0) << " " << ret(i + 1, 1) <<"\n";
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

  for (int32 i = 0; i != ret_edges; i++) {
    const int32 parent = ret(i, 0), child = ret(i, 1);
    assert(parent < (start_tip + start_node + 1));
    assert(child < (start_tip + start_node + 1));
    assert(child > 0);
    assert(new_no[parent]);
    if (!new_no[child]) {
      new_no[child] = next_no++;
    }


    ret(i, 0) = new_no[parent];
    ret(i, 1) = new_no[child];
  }


  return ret(Range(0, ret_edges - 1), _);

}
