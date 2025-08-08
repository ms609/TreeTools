#include <Rcpp/Lightest>
using namespace Rcpp;

#include "../inst/include/TreeTools/assert.h" /* for ASSERT */
#include "../inst/include/TreeTools/SplitList.h"
#include "../inst/include/TreeTools/renumber_tree.h"
using namespace TreeTools;

inline void insert_ancestor(const int16 tip, const int16 *next_node,
                            std::array<int16, SL_MAX_TIPS + SL_MAX_SPLITS>& parent,
                            std::array<int16, SL_MAX_TIPS>& patriarch) {
  if (patriarch[tip]) {
    parent[patriarch[tip]] = *next_node;
  } else {
    parent[tip] = *next_node;
  }
  patriarch[tip] = *next_node;
}

// [[Rcpp::export]]
IntegerMatrix splits_to_edge(const RawMatrix splits, const IntegerVector nTip) {
  if (double(nTip[0]) > double(std::numeric_limits<int16>::max())) {
    Rcpp::stop("This many tips are not (yet) supported.");
  }
  const int16 n_tip = int16(nTip[0]);
  if (splits.nrow() == 0) {
    IntegerMatrix ret(n_tip, 2);
    for (int i = 0; i < n_tip; ++i) {
      ret(i, 0) = n_tip + 1;
      ret(i, 1) = i + 1;
    }
    return ret;
  }
  const SplitList x(splits);
  alignas(64) std::array<int16, SL_MAX_TIPS + SL_MAX_SPLITS> parent{};
  alignas(64) std::array<int16, SL_MAX_TIPS> patriarch{};

  std::array<int16, SL_MAX_SPLITS> split_order;
  std::iota(split_order.begin(), split_order.begin() + x.n_splits, 0);
  std::sort(split_order.begin(), split_order.begin() + x.n_splits,
            [&in_split = x.in_split](int16 a, int16 b) {
              return in_split[a] > in_split[b];
            });

  int16 next_node = n_tip;
  for (int16 split = x.n_splits; split--; ) { // Work back through splits
    if (split > 0) {
      __builtin_prefetch(&x.state[split_order[split - 1]][0], 0, 3);
    }
    
    for (int16 bin = 0; bin < x.n_bins; ++bin) {
      splitbit chunk = x.state[split_order[split]][bin];
      if (!chunk) continue;
      const int16 base_tip = bin * SL_BIN_SIZE;
      while (chunk) {
        const int16 bin_tip = __builtin_ctzll(chunk); // count trailing zeros
        const int16 tip = base_tip + bin_tip;
        insert_ancestor(tip, &next_node, parent, patriarch);
        chunk &= chunk - 1; // clear lowest set bit
      }
    }
    ++next_node;
  }
  for (int16 tip = 0; tip < n_tip; ++tip) {
    insert_ancestor(tip, &next_node, parent, patriarch);
  }

  const int16 n_edge = n_tip + x.n_splits;
  IntegerVector edge1(n_edge), edge2(n_edge);
  for (int16 i = 0; i < n_edge; ++i) {
    edge1[i] = parent[i] + 1;
    edge2[i] = i + 1;
  }

  return preorder_edges_and_nodes(edge1, edge2);
}
