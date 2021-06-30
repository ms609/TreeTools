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
IntegerMatrix drop_tip (const IntegerMatrix preorder, const IntegerVector drop) {

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
