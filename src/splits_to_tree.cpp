#include <Rcpp.h>
using namespace Rcpp;

#include "../inst/include/TreeTools/SplitList.h"
#include "../inst/include/TreeTools/renumber_tree.h"
using namespace TreeTools;

inline void swap(int16 *a, int16 *b) {
  const int16 temp = *a;
  *a = *b;
  *b = temp;
}

inline void insertion_sort_by_smallest(int16* arr, const int16 arr_len,
                                       const int16* sort_by) {
  assert(arr_len > 0);
  switch (arr_len) {
  // case 0: return;
  case 1: return;
  case 2:
    if (sort_by[arr[0]] > sort_by[arr[1]]) {
      swap(&arr[0], &arr[1]);
    }
    return;
  }

  for (int16 i = 1; i != arr_len; ++i) {
    const int16
      tmp = arr[i],
      key = sort_by[tmp]
    ;
    int16 j = i;
    while (j && sort_by[arr[j - 1]] > key) {
      arr[j] = arr[j - 1];
      --j;
    }
    arr[j] = tmp;
  }
}

inline void insert_ancestor(const int16 tip, const int16 next_node,
                            int16* parent) {
  if (parent[tip] && parent[tip] != next_node) {
    insert_ancestor(parent[tip], next_node, parent);
  } else {
    parent[tip] = next_node;
    Rcout << "Parent of " << tip << " set to " << next_node << ".\n";
  }
}

// [[Rcpp::export]]
IntegerMatrix splits_to_edge(const RawMatrix splits, const IntegerVector nTip) {
  const SplitList x(splits);
  const int16
    n_tip = nTip[0],
    n_edge = n_tip + x.n_splits + 1
  ;
  int16 parent[SL_MAX_TIPS + SL_MAX_SPLITS]{};

  int16 split_order[SL_MAX_SPLITS];
  for (int16 i = x.n_splits; i--; ) {
    split_order[i] = i;
  }
  Rcout << "\n\nsplits_to_edge: " << x.n_splits << " splits loaded.\n";
  insertion_sort_by_smallest(split_order, x.n_splits, x.in_split);

  int16 next_node = n_tip;
  for (int16 split = x.n_splits; split--; ) {
    for (int16 bin = x.n_bins; --bin; ) {
      const splitbit chunk = x.state[split][bin];
      for (int16 bin_tip = SL_BIN_SIZE; bin_tip--; ) {
        const int16 tip = bin_tip + (bin * SL_BIN_SIZE);
        if (chunk & powers_of_two[tip]) {
          Rcout << "  Split " << split << " contains " << tip << ".\n";
          insert_ancestor(tip, next_node, parent);
        }
      }
    }
    ++next_node;
  }
  for (int16 tip = n_tip; tip--; ) {
    insert_ancestor(tip, next_node, parent);
  }

  IntegerVector edge1(n_edge), edge2(n_edge);
  for (int16 i = n_edge; i--; ) {
    edge1[i] = parent[i] + 1;
    edge2[i] = i + 1;
  }

  return preorder_edges_and_nodes(edge1, edge2);
}
