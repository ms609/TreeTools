#include <Rcpp.h>
using namespace Rcpp;

#include "../inst/include/TreeTools/ClusterTable.h" /* for ClusterTable */
#include "../inst/include/TreeTools/information.h" /* for split_XXX_info */
using namespace TreeTools;

#include <array> /* for array */
#include <vector> /* for vector */

// trees is a list of objects of class phylo, all with the same tip labels
// (try RenumberTips(trees, trees[[1]]))
// [[Rcpp::export]]
LogicalMatrix consensus_tree(const List trees, const NumericVector p) {
  int16
    v = 0, w = 0,
    L, R, N, W,
    L_j, R_j, N_j, W_j
  ;
  const int16
    n_trees = trees.length(),
    frac_thresh = int16(n_trees * p[0]) + 1,
    thresh = frac_thresh > n_trees ? n_trees : frac_thresh
  ;

  std::vector<ClusterTable> tables;
  tables.reserve(n_trees);
  for (int16 i = n_trees; i--; ) {
    tables.emplace_back(ClusterTable(List(trees(i))));
  }

  const int16 n_tip = tables[0].N();

  std::array<int16, CT_STACK_SIZE * CT_MAX_LEAVES> S;
  std::array<int16, CT_MAX_LEAVES> split_count;
  std::array<int16, CT_MAX_LEAVES> final_l;
  std::array<int16, CT_MAX_LEAVES> final_r;
  std::array<int16, CT_MAX_LEAVES> final_tree;

  const std::size_t ntip_3 = n_tip - 3;
  int16
    i = 0,
    splits_found = 0
  ;
  do {
    if (tables[i].NOSWX(ntip_3)) {
      continue;
    }

    std::vector<int16> split_l(n_tip);
    std::vector<int16> split_r(n_tip);
    std::fill(split_count.begin(), split_count.begin() + n_tip, 1);

    for (int16 j = i + 1; j != n_trees; j++) {

      tables[i].CLEAR();

      tables[j].TRESET();
      tables[j].READT(&v, &w);

      int16 j_pos = 0, Spos = 0; // Empty the stack S

      do {
        if (CT_IS_LEAF(v)) {
          CT_PUSH(tables[i].ENCODE(v), tables[i].ENCODE(v), 1, 1);
        } else {
          CT_POP(L, R, N, W_j);
          W = 1 + W_j;
          w = w - W_j;
          while (w) {
            CT_POP(L_j, R_j, N_j, W_j);
            if (L_j < L) L = L_j;
            if (R_j > R) R = R_j;
            N = N + N_j;
            W = W + W_j;
            w = w - W_j;
          };
          CT_PUSH(L, R, N, W);

          ++j_pos;
          if (tables[j].GETSWX(&j_pos)) {
            // Split has already been counted; next!
          } else {
            if (N == R - L + 1) { // L..R is contiguous, and must be tested
              if (tables[i].CLUSTONL(&L, &R)) {
                tables[j].SETSWX(&j_pos);
                assert(L > 0);
                ++split_count[L - 1];
                split_l[L - 1] = L;
                split_r[L - 1] = R;
                assert(split_size[L - 1] > 0);
              } else if (tables[i].CLUSTONR(&L, &R)) {
                tables[j].SETSWX(&j_pos);
                assert(R > 0);
                ++split_count[R - 1];
                split_l[R - 1] = L;
                split_r[R - 1] = R;
              }
            }
          }
        }
        tables[j].NVERTEX_short(&v, &w);
      } while (v);
    }

    for (int16 k = n_tip; k--; ) {
      if (split_count[k] >= thresh) {
        final_tree[splits_found] = i;
        final_l[splits_found] = split_l[k];
        final_r[splits_found] = split_r[k];
        ++splits_found;
        // If we have a perfectly resolved tree, break.
        if (splits_found == n_tip - 3) {
          goto rtn;
        }
      }
    }
  } while (i++ != n_trees - thresh); // All clades in p% consensus must occur in first q% of trees.
  rtn:
  LogicalMatrix ret(splits_found, n_tip);
  for (int16 i = splits_found; i--; ) {
    ClusterTable tree = tables[final_tree[i]];
    for (int16 j = final_l[i]; j != final_r[i] + 1; ++j) {
      ret(i, tree.ENCODE(j) - 1) = true;
    }
  }
  return ret;
}

// trees is a list of objects of class phylo, all with the same tip labels
// (try RenumberTips(trees, trees[[1]]))
// [[Rcpp::export]]
double consensus_info (const List trees, const LogicalVector phylo) {

  int16 v = 0, w = 0,
    L, R, N, W,
    L_j, R_j, N_j, W_j
  ;
  const int16 n_trees = trees.length();

  std::vector<ClusterTable> tables;
  tables.reserve(n_trees);
  for (int16 i = n_trees; i--; ) {
    tables.emplace_back(ClusterTable(List(trees(i))));
  }

  const int16
    n_tip = tables[0].N(),
      thresh = (n_trees / 2) + 1
  ;

  const bool phylo_info = phylo[0];

  std::array<int16, CT_STACK_SIZE * CT_MAX_LEAVES> S;
  std::array<int16, CT_MAX_LEAVES> split_count;

  double info = 0;

  const std::size_t ntip_3 = n_tip - 3;
  // All clades in 50% consensus must occur in first 50% of trees.
  for (int16 i = 0; i != thresh; i++) {
    if (tables[i].NOSWX(ntip_3)) {
      continue;
    }

    std::vector<int16> split_size(n_tip);
    std::fill(split_count.begin(), split_count.begin() + n_tip, 1);

    for (int16 j = i + 1; j != n_trees; j++) {

      tables[i].CLEAR();

      tables[j].TRESET();
      tables[j].READT(&v, &w);

      int16 j_pos = 0, Spos = 0; // Empty the stack S

      do {
        if (CT_IS_LEAF(v)) {
          CT_PUSH(tables[i].ENCODE(v), tables[i].ENCODE(v), 1, 1);
        } else {
          CT_POP(L, R, N, W_j);
          W = 1 + W_j;
          w = w - W_j;
          while (w) {
            CT_POP(L_j, R_j, N_j, W_j);
            if (L_j < L) L = L_j;
            if (R_j > R) R = R_j;
            N = N + N_j;
            W = W + W_j;
            w = w - W_j;
          };
          CT_PUSH(L, R, N, W);

          ++j_pos;
          if (tables[j].GETSWX(&j_pos)) {
            // Split has already been counted; next!
          } else {
            if (N == R - L + 1) { // L..R is contiguous, and must be tested
              if (tables[i].CLUSTONL(&L, &R)) {
                tables[j].SETSWX(&j_pos);
                assert(L > 0);
                ++split_count[L - 1];
                if (!split_size[L - 1]) {
                  split_size[L - 1] = N;
                }
                assert(split_size[L - 1] > 0);
              } else if (tables[i].CLUSTONR(&L, &R)) {
                tables[j].SETSWX(&j_pos);
                assert(R > 0);
                ++split_count[R - 1];
                if (!split_size[R - 1]) {
                  split_size[R - 1] = N;
                }
                assert(split_size[R - 1] > 0);
              }
            }
          }
        }
        tables[j].NVERTEX_short(&v, &w);
      } while (v);
    }

    int16 splits_found = 0;
    for (int16 k = n_tip; k--; ) {
      if (split_count[k] >= thresh) {
        ++splits_found;
        if (phylo_info) {
          info += split_phylo_info(split_size[k], &n_tip,
                                   split_count[k] / double(n_trees));
        } else {
          info += split_clust_info(split_size[k], &n_tip,
                                   split_count[k] / double(n_trees));
        }
        // If we have a perfectly resolved tree, break.
        if (splits_found == n_tip - 3) {
          return phylo_info ? info : info * n_tip;
        }
      }
    }
  }

  // Convert clustering entropy to *total* information
  return phylo_info ? info : info * n_tip;
}
