#include <Rcpp/Lightest>
#include <Rcpp/Interrupt.h> /* for Rcpp::checkUserInterrupt() */
using namespace Rcpp;

#include "../inst/include/TreeTools/assert.h" /* for ASSERT */
#include "../inst/include/TreeTools/ClusterTable.h" /* for ClusterTable */

#include <algorithm> /* for fill */
#include <array> /* for array */
#include <chrono> /* for steady_clock (interrupt timing) */
#include <cstdint> /* for uint64_t (split hashing) */
#include <string> /* for string (hash key) */
#include <unordered_map> /* for unordered_map */
#include <vector> /* for vector */

using TreeTools::ct_stack_threshold;
using TreeTools::ct_max_leaves_heap;
using TreeTools::ClusterTable;

struct StackEntry { int32 L, R, N, W; };

// A distinct split's "witness": tree index t and the [L, R] encode-range that
// the split spans in THAT tree's own ClusterTable encoding (where every cluster
// is contiguous). The split's leaf set is reconstructed on demand as
// {tables[t].DECODE(j) : L <= j <= R}. Lets the hashed counter defer building
// the packed bit pattern until thresholding, so non-surviving splits (the bulk,
// on low-concordance input) are never materialised.
struct SplitWitness { int32 t, L, R; };

// Reconstruct a split's packed bit pattern from its witness into row `dst`
// (nbin bytes, pre-zeroed by caller). tables is non-const (DECODE mutates none
// of the observable state but is not const-qualified).
inline void materialise_witness(std::vector<ClusterTable>& tables,
                                const SplitWitness& w, Rbyte* const dst) {
  ClusterTable& tab = tables[w.t];
  for (int32 j = w.L; j <= w.R; ++j) {
    const int32 leaf_idx = tab.DECODE(j) - 1; // 0-based original leaf id
    dst[leaf_idx >> 3] |= static_cast<Rbyte>(1 << (leaf_idx & 7));
  }
}

// Throttled (~1 s) user-interrupt check.
inline void throttled_interrupt(
    std::chrono::steady_clock::time_point& last) {
  const auto now = std::chrono::steady_clock::now();
  if (std::chrono::duration_cast<std::chrono::seconds>(now - last).count()
        >= 1) {
    last = now;
    Rcpp::checkUserInterrupt();
  }
}

// ---------------------------------------------------------------------------
// Shared postorder cluster-enumeration primitive.
//
// Walks `tree` in (short) postorder, computing each internal node's cluster
// <L, R> (leaf count N, weight W) under the leaf encoding of `ref`, and calls
// `visit(L, R, N, j_pos)` for each internal node (j_pos = 1-based index in the
// traversal).  `NVERTEX_short` skips the three trivial top vertices (root,
// ingroup, sentinel), so only non-trivial clusters are visited.  This single
// implementation backs the consensus count, the exact split-frequency count,
// and (with ref == tree) the all-splits count, so the delicate stack
// arithmetic lives in exactly one place.
template <typename Visit>
inline void for_each_internal_node(ClusterTable& ref, ClusterTable& tree,
                                   StackEntry* const S_start, Visit&& visit) {
  int32 v = 0, w = 0, L, R, N, W;
  tree.TRESET();
  tree.READT(&v, &w);
  int32 j_pos = 0;
  StackEntry* S_top = S_start; // Empty the stack S
  do {
    if (tree.is_leaf(v)) {
      const auto enc_v = ref.ENCODE(v);
      *S_top++ = {enc_v, enc_v, 1, 1};
    } else {
      const StackEntry& entry = *--S_top;
      L = entry.L; R = entry.R; N = entry.N;
      W = 1 + entry.W;
      w -= entry.W;
      while (w) {
        const StackEntry& next = *--S_top;
        L = std::min(L, next.L); // Faster than ternary operator
        R = std::max(R, next.R);
        N += next.N;
        W += next.W;
        w -= next.W;
      }
      *S_top++ = {L, R, N, W};
      ++j_pos;
      visit(L, R, N, j_pos);
    }
    tree.NVERTEX_short(&v, &w);
  } while (v);
}

// ---------------------------------------------------------------------------
// Exact single-pass count of every distinct non-trivial split.
//
// One O(k * sum-of-cluster-sizes) = O(k n h) pass: each tree's clusters are
// enumerated under its own encoding (so each is a contiguous L..R), the leaf
// set is packed into a bit pattern, and that exact pattern is the hash-map key.
// Counts accumulate directly (each tree contributes 1 per cluster it holds).
// Deterministic; no collision risk.  Backs both `split_frequencies(exact=TRUE)`
// and the majority/threshold `Consensus()`.
template <typename StackContainer>
void count_splits_exact(std::vector<ClusterTable>& tables, const int32 n_tip,
                        const int32 nbin, StackContainer& S,
                        std::vector<std::vector<Rbyte>>& split_patterns,
                        std::vector<int32>& counts) {
  const int32 n_trees = int32(tables.size());
  const int32 ntip_3 = n_tip - 3;
  std::unordered_map<std::string, int32> split_map;
  split_map.reserve((ntip_3 > 0 ? ntip_3 : 1) * 2);
  std::string key(nbin, '\0');
  StackEntry* const S_start = S.data();
  auto lastInterrupt = std::chrono::steady_clock::now();

  for (int32 t = 0; t < n_trees; ++t) {
    throttled_interrupt(lastInterrupt);
    ClusterTable& tree = tables[t];
    for_each_internal_node(tree, tree, S_start,
      [&tree, &key, &split_map, &split_patterns, &counts]
      (int32 L, int32 R, int32 /* N */, int32 /* j_pos */) {
        std::fill(key.begin(), key.end(), '\0');
        for (int32 j = L; j <= R; ++j) {
          const int32 leaf_idx = tree.DECODE(j) - 1; // 0-based
          key[leaf_idx >> 3] |= static_cast<char>(1 << (leaf_idx & 7));
        }
        auto it = split_map.find(key);
        if (it == split_map.end()) {
          split_map.emplace(key, int32(split_patterns.size()));
          split_patterns.emplace_back(key.begin(), key.end());
          counts.push_back(1);
        } else {
          ++counts[it->second];
        }
      });
  }
}

// Forward declaration: the O(kn) hashed counter (defined below) is the default
// for the majority/threshold consensus.  Non-dependent name in the template
// below, so it must be declared first.
void count_splits_hashed(std::vector<ClusterTable>& tables, const int32 n_tip,
                         std::vector<SplitWitness>& witnesses,
                         std::vector<int32>& counts);

// ---------------------------------------------------------------------------
// Consensus tree.
//
// Strict (p = 1, thresh == n_trees) keeps its already-optimal single-reference
// path over the first tree.  Majority / threshold (0.5 <= p < 1) counts every
// split's frequency in one pass and keeps those reaching the threshold: any two
// such splits each occur in > k/2 trees, so they co-occur in some tree and are
// pairwise (hence globally) compatible, forming a valid tree directly.  The
// count is hashed (O(kn), probabilistic) by default, or exact (deterministic,
// O(k.n.height)) when `exact` is set.
template<typename StackContainer>
RawMatrix calc_consensus_tree(
  const List& trees,
  const NumericVector& p,
  const bool exact,
  StackContainer& S
) {
  const int32 n_trees = trees.length();
  const int32 frac_thresh = int32(n_trees * p[0]) + 1;
  const int32 thresh = frac_thresh > n_trees ? n_trees : frac_thresh;

  std::vector<ClusterTable> tables;
  tables.reserve(n_trees);
  for (int32 i = 0; i < n_trees; ++i) {
    tables.emplace_back(ClusterTable(Rcpp::List(trees(i))));
  }

  const int32 n_tip   = tables[0].N();
  const int32 ntip_3  = n_tip - 3;
  const int32 nbin    = (n_tip + 7) / 8;  // bytes per row in packed output

  StackEntry *const S_start = S.data();
  RawMatrix ret(ntip_3, nbin);
  int32 splits_found = 0;

  if (thresh >= n_trees) {
    // ---- Strict consensus: single reference (tree 0) ----------------------
    int32* split_count;
    std::array<int32, ct_stack_threshold> split_stack;
    std::vector<int32> split_heap;
    if (n_tip <= ct_stack_threshold) {
      split_count = split_stack.data();
    } else {
      split_heap.resize(n_tip);
      split_count = split_heap.data();
    }
    std::fill(split_count, split_count + n_tip, 1); // tree 0 holds its clusters

    auto lastInterrupt = std::chrono::steady_clock::now();
    for (int32 j = 1; j < n_trees; ++j) {
      throttled_interrupt(lastInterrupt);
      ASSERT(tables[j].N() == n_tip);
      for_each_internal_node(tables[0], tables[j], S_start,
        [&tables, &split_count](int32 cl_L, int32 cl_R, int32 cl_N, int32) {
          if (cl_N == cl_R - cl_L + 1) { // contiguous: testable against ref
            if (tables[0].CLUSTONL(cl_L, cl_R)) {
              ASSERT(cl_L > 0);
              ++split_count[cl_L - 1];
            } else if (tables[0].CLUSTONR(cl_L, cl_R)) {
              ASSERT(cl_R > 0);
              ++split_count[cl_R - 1];
            }
          }
        });
    }
    // Pack reference clusters present in every tree.
    for (int32 k = 0; k < n_tip; ++k) {
      if (split_count[k] >= thresh) {
        const int32 start = tables[0].X_left(k + 1);
        const int32 end   = tables[0].X_right(k + 1);
        if (start == 0 && end == 0) continue; // no cluster at this row
        for (int32 j = start; j <= end; ++j) {
          const int32 leaf_idx = tables[0].DECODE(j) - 1;
          const int32 byte_idx = leaf_idx >> 3;
          const int32 bit_idx  = leaf_idx & 7;
          Rbyte* col_ptr = &ret(0, byte_idx);
          col_ptr[splits_found] |= (Rbyte(1) << bit_idx);
        }
        ++splits_found;
        if (splits_found == ntip_3) return ret;
      }
    }
  } else if (exact) {
    // ---- Majority / threshold, deterministic: count then threshold --------
    std::vector<std::vector<Rbyte>> split_patterns;
    std::vector<int32> counts;
    count_splits_exact(tables, n_tip, nbin, S, split_patterns, counts);
    const int32 n_distinct = int32(split_patterns.size());
    for (int32 i = 0; i < n_distinct; ++i) {
      if (counts[i] >= thresh) {
        for (int32 c = 0; c < nbin; ++c) {
          ret(splits_found, c) = split_patterns[i][c];
        }
        ++splits_found;
        if (splits_found == ntip_3) return ret;
      }
    }
  } else {
    // ---- Majority / threshold, hashed: count, then materialise ONLY the
    //      splits that reach the threshold (the rest are never built). --------
    std::vector<SplitWitness> witnesses;
    std::vector<int32> counts;
    count_splits_hashed(tables, n_tip, witnesses, counts);
    const int32 n_distinct = int32(witnesses.size());
    std::vector<Rbyte> row(nbin);
    for (int32 i = 0; i < n_distinct; ++i) {
      if (counts[i] >= thresh) {
        std::fill(row.begin(), row.end(), Rbyte(0));
        materialise_witness(tables, witnesses[i], row.data());
        for (int32 c = 0; c < nbin; ++c) {
          ret(splits_found, c) = row[c];
        }
        ++splits_found;
        if (splits_found == ntip_3) return ret;
      }
    }
  }

  return (splits_found == 0) ? RawMatrix(0, nbin) :
         (splits_found < ntip_3) ? ret(Range(0, splits_found - 1), _) : ret;
}

// ---------------------------------------------------------------------------
// Hashed split frequencies (the fast default for split_frequencies).
//
// A single O(kn) pass: each non-trivial cluster is identified by a 128-bit
// subtree hash = the (order-independent) sum of its leaves' fixed splitmix64
// hashes, so the same split in different trees hashes identically without an
// O(cluster size) key.  Counts accumulate directly; materialisation of the
// packed bit pattern is DEFERRED — each distinct split keeps only a witness
// (tree + encode-range), so the bulk of splits that never reach the consensus
// threshold are never materialised (the dominant cost at high n).  Exactness is
// probabilistic (a 128-bit collision, ~1e-30, would conflate two splits).
inline uint64_t splitmix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

struct HashEntry { int32 L, R, W; uint64_t lo, hi; };
struct Hash128 {
  uint64_t lo, hi;
  bool operator==(const Hash128& o) const noexcept {
    return lo == o.lo && hi == o.hi;
  }
};
struct Hash128Hasher {
  std::size_t operator()(const Hash128& k) const noexcept {
    return std::size_t(k.lo ^ (k.hi * 0x9e3779b97f4a7c15ULL));
  }
};

void count_splits_hashed(std::vector<ClusterTable>& tables, const int32 n_tip,
                         std::vector<SplitWitness>& witnesses,
                         std::vector<int32>& counts) {
  const int32 n_trees = int32(tables.size());
  const int32 ntip_3 = n_tip - 3;

  // Fixed per-leaf 128-bit hashes, keyed by original (1-based) leaf id.
  std::vector<uint64_t> leaf_lo(n_tip + 1), leaf_hi(n_tip + 1);
  for (int32 i = 1; i <= n_tip; ++i) {
    leaf_lo[i] = splitmix64(uint64_t(i));
    leaf_hi[i] = splitmix64(uint64_t(i) + 0x9e3779b97f4a7c15ULL);
  }

  HashEntry* S_start;
  std::array<HashEntry, ct_stack_threshold> hash_stack;
  std::vector<HashEntry> hash_heap;
  if (n_tip <= ct_stack_threshold) {
    S_start = hash_stack.data();
  } else {
    hash_heap.resize(n_tip);
    S_start = hash_heap.data();
  }

  std::unordered_map<Hash128, int32, Hash128Hasher> split_map;
  split_map.reserve((ntip_3 > 0 ? ntip_3 : 1) * 2);

  auto lastInterrupt = std::chrono::steady_clock::now();

  for (int32 t = 0; t < n_trees; ++t) {
    throttled_interrupt(lastInterrupt);
    ClusterTable& tree = tables[t];
    int32 v = 0, w = 0, L, R, W;
    tree.TRESET();
    tree.READT(&v, &w);
    HashEntry* S_top = S_start;
    do {
      if (tree.is_leaf(v)) {
        const int32 enc_v = tree.ENCODE(v);
        *S_top++ = {enc_v, enc_v, 1, leaf_lo[v], leaf_hi[v]};
      } else {
        const HashEntry& entry = *--S_top;
        L = entry.L; R = entry.R;
        W = 1 + entry.W;
        uint64_t hlo = entry.lo, hhi = entry.hi;
        w -= entry.W;
        while (w) {
          const HashEntry& next = *--S_top;
          L = std::min(L, next.L);
          R = std::max(R, next.R);
          W += next.W;
          hlo += next.lo;
          hhi += next.hi;
          w -= next.W;
        }
        *S_top++ = {L, R, W, hlo, hhi};

        const Hash128 hkey{hlo, hhi};
        auto it = split_map.find(hkey);
        if (it == split_map.end()) {
          split_map.emplace(hkey, int32(witnesses.size()));
          // Defer materialisation: record where this split lives (tree t,
          // encode-range [L, R]) and rebuild the pattern later only if needed.
          witnesses.push_back({t, L, R});
          counts.push_back(1);
        } else {
          ++counts[it->second];
        }
      }
      tree.NVERTEX_short(&v, &w);
    } while (v);
  }
}

// Assemble a split-frequency result List from collected patterns + counts.
inline List frequencies_list(
    const std::vector<std::vector<Rbyte>>& split_patterns,
    const std::vector<int32>& counts, const int32 nbin) {
  const int32 splits_found = int32(split_patterns.size());
  RawMatrix ret(splits_found, nbin);
  for (int32 r = 0; r < splits_found; ++r) {
    for (int32 c = 0; c < nbin; ++c) {
      ret(r, c) = split_patterns[r][c];
    }
  }
  IntegerVector count_vec(counts.begin(), counts.end());
  return List::create(Named("splits") = ret, Named("counts") = count_vec);
}

// As frequencies_list, but materialising each split's pattern from its witness
// (split_frequencies needs every distinct split, so all are built — same total
// work as eager, just deferred, with less peak memory: no vector-of-patterns).
inline List frequencies_list_from_witnesses(
    std::vector<ClusterTable>& tables,
    const std::vector<SplitWitness>& witnesses,
    const std::vector<int32>& counts, const int32 nbin) {
  const int32 splits_found = int32(witnesses.size());
  RawMatrix ret(splits_found, nbin);
  std::vector<Rbyte> row(nbin);
  for (int32 r = 0; r < splits_found; ++r) {
    std::fill(row.begin(), row.end(), Rbyte(0));
    materialise_witness(tables, witnesses[r], row.data());
    for (int32 c = 0; c < nbin; ++c) {
      ret(r, c) = row[c];
    }
  }
  IntegerVector count_vec(counts.begin(), counts.end());
  return List::create(Named("splits") = ret, Named("counts") = count_vec);
}

List calc_split_frequencies_hashed(const List& trees, const int32 n_tip) {
  const int32 n_trees = trees.length();
  std::vector<ClusterTable> tables;
  tables.reserve(n_trees);
  for (int32 i = 0; i < n_trees; ++i) {
    tables.emplace_back(ClusterTable(Rcpp::List(trees(i))));
  }
  const int32 nbin = (n_tip + 7) / 8;
  std::vector<SplitWitness> witnesses;
  std::vector<int32> counts;
  count_splits_hashed(tables, n_tip, witnesses, counts);
  return frequencies_list_from_witnesses(tables, witnesses, counts, nbin);
}

// Exact split frequencies: the same single-pass count, returned in full.
template<typename StackContainer>
List calc_split_frequencies_exact(const List& trees, StackContainer& S) {
  const int32 n_trees = trees.length();
  std::vector<ClusterTable> tables;
  tables.reserve(n_trees);
  for (int32 i = 0; i < n_trees; ++i) {
    tables.emplace_back(ClusterTable(Rcpp::List(trees(i))));
  }
  const int32 n_tip = tables[0].N();
  const int32 nbin  = (n_tip + 7) / 8;

  std::vector<std::vector<Rbyte>> split_patterns;
  std::vector<int32> counts;
  count_splits_exact(tables, n_tip, nbin, S, split_patterns, counts);
  return frequencies_list(split_patterns, counts, nbin);
}

// ---------------------------------------------------------------------------
// Exports

// [[Rcpp::export]]
List split_frequencies(const List trees, const bool exact = false) {
  try {
    ClusterTable temp_table(Rcpp::List(trees(0)));
    const int32 n_tip = temp_table.N();

    if (!exact) {
      return calc_split_frequencies_hashed(trees, n_tip);
    }
    if (n_tip <= ct_stack_threshold) {
      std::array<StackEntry, ct_stack_threshold> S;
      return calc_split_frequencies_exact(trees, S);
    } else {
      std::vector<StackEntry> S(n_tip);
      return calc_split_frequencies_exact(trees, S);
    }
  } catch(const std::exception& e) {
    Rcpp::stop(e.what());
  }

  ASSERT(false && "Unreachable code in split_frequencies");
  return List();
}

// trees is a list of objects of class phylo, all with the same tip labels
// (try RenumberTips(trees, trees[[1]]))
// [[Rcpp::export]]
RawMatrix consensus_tree(const List trees, const NumericVector p,
                         const bool exact = false) {
  try {
    ClusterTable temp_table(Rcpp::List(trees(0)));
    const int32 n_tip = temp_table.N();

    if (n_tip <= ct_stack_threshold) {
      std::array<StackEntry, ct_stack_threshold> S;
      return calc_consensus_tree(trees, p, exact, S);
    } else {
      std::vector<StackEntry> S(n_tip);
      return calc_consensus_tree(trees, p, exact, S);
    }
  } catch(const std::exception& e) {
    Rcpp::stop(e.what());
  }

  ASSERT(false && "Unreachable code in consensus_tree");
  return RawMatrix(0, 0);
}
