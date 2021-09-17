#include <Rcpp.h>
using namespace Rcpp;

#include "../inst/include/TreeTools/SplitList.h"
#include "../inst/include/TreeTools/renumber_tree.h"
using namespace TreeTools;

TREETOOLS_SPLITLIST_INIT

inline void insertion_sort_by_largest(int16* arr, const int16 arr_len,
                                       const int16* sort_by) {
  assert(arr_len > 0);
  switch (arr_len) {
  // case 0: return;
  case 1: return;
  case 2:
    if (sort_by[arr[0]] < sort_by[arr[1]]) {
      const int16 tmp = arr[0];
      arr[0] = arr[1];
      arr[1] = tmp;
    }
    return;
  }

  for (int16 i = 1; i != arr_len; ++i) {
    const int16
      tmp = arr[i],
      key = sort_by[tmp]
    ;
    int16 j = i;
    while (j && sort_by[arr[j - 1]] < key) {
      arr[j] = arr[j - 1];
      --j;
    }
    arr[j] = tmp;
  }
}

inline void insert_ancestor(const int16 tip, const int16 *next_node,
                            int16 (&parent)[SL_MAX_TIPS + SL_MAX_SPLITS],
                            int16 (&patriarch)[SL_MAX_TIPS]) {
  if (patriarch[tip]) {
    parent[patriarch[tip]] = *next_node;
  } else {
    parent[tip] = *next_node;
  }
  patriarch[tip] = *next_node;
}

// [[Rcpp::export]]
IntegerMatrix splits_to_edge(const RawMatrix splits, const IntegerVector nTip) {
  const SplitList x(splits);
  const int16 n_tip = nTip[0];
  int16 parent[SL_MAX_TIPS + SL_MAX_SPLITS]{};
  int16 patriarch[SL_MAX_TIPS]{};

  int16 split_order[SL_MAX_SPLITS];
  for (int16 i = x.n_splits; i--; ) {
    split_order[i] = i;
  }
  // Rcout << "\n\nsplits_to_edge: " << x.n_splits << " splits loaded.\n";
  insertion_sort_by_largest(split_order, x.n_splits, x.in_split);

  int16 next_node = n_tip;
  for (int16 split = x.n_splits; split--; ) {
    for (int16 bin = x.n_bins; bin--; ) {
      const splitbit chunk = x.state[split_order[split]][bin];
      for (int16 bin_tip = SL_BIN_SIZE; bin_tip--; ) {
        const int16 tip = bin_tip + (bin * SL_BIN_SIZE);
        if (chunk & powers_of_two[bin_tip]) {
          insert_ancestor(tip, &next_node, parent, patriarch);
        }
      }
    }
    ++next_node;
  }
  for (int16 tip = n_tip; tip--; ) {
    insert_ancestor(tip, &next_node, parent, patriarch);
  }

  const int16 n_edge = n_tip + x.n_splits;
  IntegerVector edge1(n_edge), edge2(n_edge);
  for (int16 i = n_edge; i--; ) {
    edge1[i] = parent[i] + 1;
    edge2[i] = i + 1;
    // Rcout << "   Edge " << i << ": " << edge1[i] << " - " << edge2[i] << ".\n";
  }

  return preorder_edges_and_nodes(edge1, edge2);
}
