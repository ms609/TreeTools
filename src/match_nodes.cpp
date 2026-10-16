#include <Rcpp/Lightest>
#include <algorithm> /* for max */
#include <climits> /* for INT_MAX */
#include <numeric> /* for iota */
#include <vector>

using namespace Rcpp;

// Match the clades of tree `x` against those of tree `table` in linear time.
//
// Table clades are encoded in a cluster table (Day 1985,
// doi:10.1007/BF01908061).  Table tips are ranked so that every clade is a
// contiguous interval [lo, hi] of ranks.  A clade of two or more tips is
// filed in row `lo` if it is the last child of its parent (or is the root),
// and in row `hi` otherwise; no two clades share a row.  A clade of `x` is
// present in `table` if and only if its tips map to `cnt == hi - lo + 1`
// distinct ranks and row `lo` or row `hi` holds [lo, hi].
//
// Nodes with a single child, and (for edge matching) table tips absent from
// `x`, mean that several table nodes can denote the same clade.  Such nodes
// form a chain, which is filed once, under the node whose incoming edge
// appears first in the edge matrix (the root, which has no edge, ranks
// last).  This reproduces the first-hit semantics of `match()` used by the
// original R implementation.
//
// Neither tree need be in any particular order: nodes are visited
// children-first using pending-child counts (Kahn's algorithm), and ranks are
// assigned from offsets within each parent, taken in edge order.
//
// @param x_edge,t_edge Two-column edge matrices of `x` and `table`.
// @param tip_match For each tip 1..n of `x`, the table node with the same
// label, or NA.  Non-NA values must be unique (enforced by the R caller).
// @param edges If `false`, return, for each node of `x` (indexed by node
// number), the matching table node; the table root is a candidate.
// If `true`, return, for each edge of `x`, the index of the matching table
// edge; table tips absent from `x` are disregarded and the root is not a
// candidate.
// @return Integer vector with NA where no match exists.
// [[Rcpp::export]]
IntegerVector match_clades(const IntegerMatrix x_edge,
                           const IntegerMatrix t_edge,
                           const IntegerVector tip_match,
                           const bool edges) {
  if (x_edge.ncol() != 2 || t_edge.ncol() != 2) {
    Rcpp::stop("Edge matrices must have two columns");
  }
  const int nx_edge = x_edge.nrow();
  const int nt_edge = t_edge.nrow();
  const int n_tip_x = tip_match.size();
  const int* x_parent = x_edge.begin();
  const int* x_child = x_parent + nx_edge;
  const int* t_parent = t_edge.begin();
  const int* t_child = t_parent + nt_edge;
  const int* tip_m = tip_match.begin();
  
  int x_nodes = n_tip_x, t_nodes = 0, x_min = 1, t_min = 1;
  for (int i = 0; i != 2 * nx_edge; ++i) {
    x_nodes = std::max(x_nodes, x_parent[i]);
    x_min = std::min(x_min, x_parent[i]);
  }
  for (int i = 0; i != 2 * nt_edge; ++i) {
    t_nodes = std::max(t_nodes, t_parent[i]);
    t_min = std::min(t_min, t_parent[i]);
  }
  if (x_min < 1 || t_min < 1) { // also catches NA
    Rcpp::stop("Node numbers must be positive integers");
  }
  const int x_size = x_nodes + 1, t_size = t_nodes + 1;
  
  // Integer scratch space, zero-initialized
  std::vector<int> block(size_t(t_size) * 9 + size_t(x_size) * 7, 0);
  int* t_par = block.data();       // parent (0: none)
  int* t_in = t_par + t_size;      // index of incoming edge
  int* t_pend = t_in + t_size;     // pending children; then child offsets
  int* t_cnt = t_pend + t_size;    // participating tips in clade
  int* t_rel = t_cnt + t_size;     // offset of clade within parent's clade
  int* t_lo = t_rel + t_size;      // first rank in clade
  int* t_ord = t_lo + t_size;      // nodes, children before parents
  int* t_head = t_ord + t_size;    // head of single-child chain
  int* t_leaf = t_head + t_size;   // 1 for tips
  int* x_par = t_leaf + t_size;
  int* x_pend = x_par + x_size;
  int* x_cnt = x_pend + x_size;    // -1: contains a tip absent from table
  int* x_lo = x_cnt + x_size;
  int* x_hi = x_lo + x_size;
  int* x_ord = x_hi + x_size;
  int* x_res = x_ord + x_size;     // matching table edge + 1 (0: none)
  
  // ---- Table: parents and child counts -----------------------------------
  for (int i = 0; i != nt_edge; ++i) {
    const int c = t_child[i], p = t_parent[i];
    if (t_par[c]) Rcpp::stop("`table` contains a node with two parents");
    t_par[c] = p;
    t_in[c] = i;
    ++t_pend[p];
  }
  
  // ---- Table: tips that take part ----------------------------------------
  // Node matching: all tips; edge matching: only those in x.
  int n_ord = 0, n_present = 0, t_root = 0, n_root = 0;
  for (int v = 1; v != t_size; ++v) {
    if (!t_par[v] && !t_pend[v]) continue; // Not in tree
    ++n_present;
    if (!t_par[v]) {
      t_root = v;
      ++n_root;
    }
    if (!t_pend[v]) {
      t_ord[n_ord++] = v;
      t_leaf[v] = 1;
      t_cnt[v] = !edges;
    }
  }
  if (n_root > 1) Rcpp::stop("`table` must be a single connected tree");
  if (edges) {
    for (int i = 0; i != n_tip_x; ++i) {
      const int m = tip_m[i];
      if (m > 0 && m < t_size && t_leaf[m]) t_cnt[m] = 1;
    }
  }
  
  // ---- Table: children-first order; clade sizes --------------------------
  for (int k = 0; k != n_ord; ++k) {
    const int v = t_ord[k];
    const int p = t_par[v];
    if (p) {
      t_cnt[p] += t_cnt[v];
      if (!--t_pend[p]) t_ord[n_ord++] = p;
    }
  }
  if (n_ord != n_present) Rcpp::stop("`table` is not a tree");
  
  // ---- Table: offset of each clade within its parent, in edge order ------
  // t_pend is now zero throughout, so can accumulate offsets.
  for (int i = 0; i != nt_edge; ++i) {
    const int c = t_child[i], p = t_parent[i];
    t_rel[c] = t_pend[p];
    t_pend[p] += t_cnt[c];
  }
  
  // ---- Table: ranks and single-child chains, parents first ---------------
  // A node continues its parent's chain if it holds all its parent's tips.
  // The chain's representative (t_best, indexed by chain head) is the member
  // whose incoming edge comes first; the root ranks last (node matching) or
  // is ineligible (edge matching).
  const int root_prio = edges ? INT_MAX : nt_edge;
  int* t_best = t_pend; // Offsets no longer needed
  for (int k = n_ord; k--; ) {
    const int v = t_ord[k];
    const int p = t_par[v];
    t_lo[v] = p ? t_lo[p] + t_rel[v] : 0;
    if (!t_cnt[v]) continue;
    const int prio = p ? t_in[v] : root_prio;
    if (p && t_cnt[v] == t_cnt[p]) {
      const int h = t_head[p];
      t_head[v] = h;
      if (prio < t_best[h]) t_best[h] = prio;
    } else {
      t_head[v] = v;
      t_best[v] = prio;
    }
  }
  
  // ---- Table: file chains in cluster table -------------------------------
  const int n_rank = t_root ? t_cnt[t_root] : 0;
  std::vector<int> rows(size_t(n_rank) * 4 + 1, INT_MAX);
  int* row_lo = rows.data();
  int* row_hi = row_lo + n_rank;
  int* row_prio = row_hi + n_rank;
  int* single = row_prio + n_rank;
  for (int k = 0; k != n_ord; ++k) {
    const int v = t_ord[k];
    const int cnt = t_cnt[v];
    if (!cnt || t_head[v] != v) continue;
    const int best = t_best[v];
    if (best == INT_MAX) continue;
    const int lo = t_lo[v];
    if (cnt == 1) {
      single[lo] = best;
    } else {
      const int p = t_par[v];
      const int hi = lo + cnt - 1;
      const int row = (!p || t_rel[v] + cnt == t_cnt[p]) ? lo : hi;
      row_lo[row] = lo;
      row_hi[row] = hi;
      row_prio[row] = best;
    }
  }
  
  // ---- x: parents and child counts ---------------------------------------
  for (int i = 0; i != nx_edge; ++i) {
    const int c = x_child[i], p = x_parent[i];
    if (x_par[c]) Rcpp::stop("`x` contains a node with two parents");
    x_par[c] = p;
    ++x_pend[p];
  }
  
  // ---- x: tips -----------------------------------------------------------
  int nx_ord = 0;
  for (int v = 1; v != x_size; ++v) {
    x_lo[v] = INT_MAX;
    x_hi[v] = -1;
    const bool present = x_par[v] || x_pend[v] || (nx_edge == 0 && v == 1);
    if (!present || x_pend[v]) continue;
    x_ord[nx_ord++] = v;
    const int m = v <= n_tip_x ? tip_m[v - 1] : NA_INTEGER;
    if (m > 0 && m < t_size && t_leaf[m] && t_cnt[m]) {
      x_lo[v] = t_lo[m];
      x_hi[v] = t_lo[m];
      x_cnt[v] = 1;
    } else {
      x_cnt[v] = -1;
    }
  }
  
  // ---- Look up clades of x, children first -------------------------------
  for (int k = 0; k != nx_ord; ++k) {
    const int v = x_ord[k];
    const int cnt = x_cnt[v];
    if (cnt > 0) {
      const int lo = x_lo[v], hi = x_hi[v];
      int prio = INT_MAX;
      if (hi - lo + 1 == cnt) {
        if (lo == hi) {
          prio = single[lo];
        } else if (row_lo[lo] == lo && row_hi[lo] == hi) {
          prio = row_prio[lo];
        } else if (row_lo[hi] == lo && row_hi[hi] == hi) {
          prio = row_prio[hi];
        }
      }
      if (prio != INT_MAX) x_res[v] = prio + 1;
    }
    const int p = x_par[v];
    if (p) {
      if (cnt < 0 || x_cnt[p] < 0) {
        x_cnt[p] = -1;
      } else {
        x_cnt[p] += cnt;
        if (x_lo[v] < x_lo[p]) x_lo[p] = x_lo[v];
        if (x_hi[v] > x_hi[p]) x_hi[p] = x_hi[v];
      }
      if (!--x_pend[p]) x_ord[nx_ord++] = p;
    }
  }
  
  if (edges) {
    IntegerVector ret(nx_edge);
    int* res = ret.begin();
    for (int i = 0; i != nx_edge; ++i) {
      const int m = x_res[x_child[i]];
      res[i] = m ? m : NA_INTEGER;
    }
    return ret;
  }
  IntegerVector ret(x_nodes);
  int* res = ret.begin();
  for (int v = 1; v != x_size; ++v) {
    const int m = x_res[v];
    res[v - 1] = m ? (m == nt_edge + 1 ? t_root : t_child[m - 1]) : NA_INTEGER;
  }
  return ret;
}

// For each internal node of a renumbered tree, the index (into the original
// `node.label` vector) of the corresponding node of the original tree.
// Both trees must have tips 1..n_tip and internal nodes numbered from
// n_tip + 1.
// [[Rcpp::export]]
IntegerVector node_label_order(const IntegerMatrix new_edge,
                               const IntegerMatrix old_edge,
                               const int n_tip) {
  IntegerVector tips(n_tip);
  std::iota(tips.begin(), tips.end(), 1);
  const IntegerVector map = match_clades(new_edge, old_edge, tips, false);
  const int n_node = map.size() - n_tip;
  if (n_node < 0) Rcpp::stop("`n_tip` exceeds number of nodes");
  IntegerVector ret(n_node);
  const int* m = map.begin() + n_tip;
  int* r = ret.begin();
  for (int i = 0; i != n_node; ++i) {
    r[i] = m[i] == NA_INTEGER ? NA_INTEGER : m[i] - n_tip;
  }
  return ret;
}
