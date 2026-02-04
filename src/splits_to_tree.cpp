#include <Rcpp/Lightest>
using namespace Rcpp;

#include "../inst/include/TreeTools/assert.h" /* for ASSERT */
#include "../inst/include/TreeTools/SplitList.h"
#include "../inst/include/TreeTools/renumber_tree.h"
using namespace TreeTools;

inline void insert_ancestor(const int32 tip, const int32 *next_node,
                            int32* parent,
                            int32* patriarch) {
  if (patriarch[tip]) {
    parent[patriarch[tip]] = *next_node;
  } else {
    parent[tip] = *next_node;
  }
  patriarch[tip] = *next_node;
}

// [[Rcpp::export]]
IntegerMatrix splits_to_edge(const RawMatrix splits, const IntegerVector nTip) {
  const int32 n_tip = int32(nTip[0]);
  if (splits.nrow() == 0) {
    IntegerMatrix ret(n_tip, 2);
    for (int i = 0; i < n_tip; ++i) {
      ret(i, 0) = n_tip + 1;
      ret(i, 1) = i + 1;
    }
    return ret;
  }
  const SplitList x(splits);
  
  // Decide whether to use stack or heap allocation based on tree size
  const bool use_heap = (n_tip > SL_MAX_TIPS) || (x.n_splits > SL_MAX_SPLITS);
  
  // Stack allocation for small trees (fast path)
  alignas(64) std::array<int32, SL_MAX_TIPS + SL_MAX_SPLITS> stack_parent{};
  alignas(64) std::array<int32, SL_MAX_TIPS> stack_patriarch{};
  
  // Heap allocation for large trees
  std::vector<int32> heap_parent;
  std::vector<int32> heap_patriarch;
  
  // Pointers to active storage
  int32* parent;
  int32* patriarch;
  
  if (use_heap) {
    const size_t parent_size = static_cast<size_t>(n_tip) + 
                               static_cast<size_t>(x.n_splits);
    heap_parent.resize(parent_size, 0);
    heap_patriarch.resize(n_tip, 0);
    parent = heap_parent.data();
    patriarch = heap_patriarch.data();
  } else {
    parent = stack_parent.data();
    patriarch = stack_patriarch.data();
  }
  
  // Allocate split_order appropriately
  std::vector<int32> split_order(x.n_splits);
  std::iota(split_order.begin(), split_order.end(), 0);
  std::sort(split_order.begin(), split_order.end(),
            [&in_split = x.in_split](int32 a, int32 b) {
              return in_split[a] > in_split[b];
            });

  int32 next_node = n_tip;
  for (int32 split = x.n_splits; split--; ) { // Work back through splits
    if (split > 0) {
      __builtin_prefetch(&x.state[split_order[split - 1]][0], 0, 3);
    }
    
    for (int32 bin = 0; bin < x.n_bins; ++bin) {
      splitbit chunk = x.state[split_order[split]][bin];
      if (!chunk) continue;
      const int32 base_tip = bin * SL_BIN_SIZE;
      while (chunk) {
        const int32 bin_tip = __builtin_ctzll(chunk); // count trailing zeros
        const int32 tip = base_tip + bin_tip;
        insert_ancestor(tip, &next_node, parent, patriarch);
        chunk &= chunk - 1; // clear lowest set bit
      }
    }
    ++next_node;
  }
  for (int32 tip = 0; tip < n_tip; ++tip) {
    insert_ancestor(tip, &next_node, parent, patriarch);
  }

  const int32 n_edge = n_tip + x.n_splits;
  IntegerVector edge1(n_edge), edge2(n_edge);
  for (int32 i = 0; i < n_edge; ++i) {
    edge1[i] = parent[i] + 1;
    edge2[i] = i + 1;
  }

  return preorder_edges_and_nodes(edge1, edge2);
}
