#include <Rcpp/Lightest>
using namespace Rcpp;

#include "../inst/include/TreeTools/assert.h" /* for ASSERT */
#include "../inst/include/TreeTools/ClusterTable.h" /* for ClusterTable */
#include "../inst/include/TreeTools/ClusterTableSIMD.h" /* for SIMD extensions */

#include <algorithm> /* for fill */
#include <array> /* for array */
#include <vector> /* for vector */
#include <chrono> /* for benchmarking */

// YOUR ORIGINAL FUNCTION - UNCHANGED
// [[Rcpp::export]]
LogicalMatrix consensus_tree(const List trees, const NumericVector p) {
  int16
  v = 0, w = 0,
    L, R, N, W,
    L_j, R_j, N_j, W_j
  ;
  const int32
  n_trees = trees.length(),
    frac_thresh = int32(n_trees * p[0]) + 1,
    thresh = frac_thresh > n_trees ? n_trees : frac_thresh
  ;
  
  std::vector<TreeTools::ClusterTable> tables;
  tables.reserve(n_trees);
  for (int32 i = n_trees; i--; ) {
    tables.emplace_back(TreeTools::ClusterTable(Rcpp::List(trees(i))));
  }
  
  const int32
  n_tip = tables[0].N(),
    ntip_3 = n_tip - 3
  ;
  
  std::array<int32, CT_STACK_SIZE * CT_MAX_LEAVES> S;
  std::array<int32, CT_MAX_LEAVES> split_count;
  
  LogicalMatrix ret(ntip_3, n_tip);
  
  int32 i = 0;
  int32 splits_found = 0;
  
  do {
    if (tables[i].NOSWX(ntip_3)) {
      continue;
    }
    
    std::fill(split_count.begin(), split_count.begin() + n_tip, 1);
    
    for (int32 j = i + 1; j != n_trees; j++) {
      
      tables[i].CLEAR();
      
      tables[j].TRESET();
      tables[j].READT(&v, &w);
      
      int16 j_pos = 0;
      int32 Spos = 0; // Empty the stack S. Used in CT_PUSH /CT_POP macros.
      
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
                ASSERT(L > 0);
                ++split_count[L - 1];
              } else if (tables[i].CLUSTONR(&L, &R)) {
                tables[j].SETSWX(&j_pos);
                ASSERT(R > 0);
                ++split_count[R - 1];
              }
            }
          }
        }
        tables[j].NVERTEX_short(&v, &w);
      } while (v);
    }
    
    for (int32 k = n_tip; k--; ) {
      if (split_count[k] >= thresh) {
        for (int32 j = tables[i].X_left(k + 1);
             j != tables[i].X_right(k + 1) + 1;
             ++j) {
          ret(splits_found, tables[i].DECODE(j) - 1) = true;
        }
        ++splits_found;
        // If we have a perfectly resolved tree, break.
        if (splits_found == ntip_3) {
          return ret;
        }
      }
    }
  } while (i++ != n_trees - thresh);
  
  return splits_found ? 
  ret(Range(0, splits_found - 1), _) :
    LogicalMatrix(0, n_tip);
}

// [[Rcpp::export]]
LogicalMatrix consensus_tree_simd(const List trees, const NumericVector p) {
  const int32
  n_trees = trees.length(),
    frac_thresh = int32(n_trees * p[0]) + 1,
    thresh = frac_thresh > n_trees ? n_trees : frac_thresh
  ;
  
  std::vector<TreeTools::ClusterTableSIMD> tables;
  tables.reserve(n_trees);
  
  for (int32 i = 0; i < n_trees; ++i) {
    const Rcpp::List& tree_i = Rcpp::List(trees(i));
    tables.emplace_back(TreeTools::ClusterTableSIMD(tree_i));
  }
  
  const int32
  n_tip = tables[0].N(),
    ntip_3 = n_tip - 3
  ;
  
  std::array<int32, CT_STACK_SIZE * CT_MAX_LEAVES> S;
  LogicalMatrix ret(ntip_3, n_tip);
  
  int32 i = 0;
  int32 splits_found = 0;
  
  int16 v = 0, w = 0, L, R, N, W, L_j, R_j, N_j, W_j;
  
  do {
    if (tables[i].NOSWX(ntip_3)) {
      continue;
    }
    
    tables[i].init_split_count_simd(n_tip);
    
    for (int32 j = i + 1; j != n_trees; j++) {
      tables[i].CLEAR();
      tables[j].TRESET();
      tables[j].READT(&v, &w);
      
      int16 j_pos = 0;
      int32 Spos = 0;
      
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
                ASSERT(L > 0);
                tables[i].increment_split_count_simd(L - 1);
              } else if (tables[i].CLUSTONR(&L, &R)) {
                tables[j].SETSWX(&j_pos);
                ASSERT(R > 0);
                tables[i].increment_split_count_simd(R - 1);
              }
            }
          }
        }
        tables[j].NVERTEX_short(&v, &w);
      } while (v);
    }
    
    tables[i].prefetch_split_data();
    // NEW: Vectorized threshold checking with early exit
    std::vector<int32> qualifying_splits;
    qualifying_splits.reserve(n_tip); // Pre-allocate
    
    tables[i].find_splits_above_threshold_simd(thresh, qualifying_splits);
    
    // Process qualifying splits
    for (int32 k : qualifying_splits) {
      int32 first_qualifying = TreeTools::simd_utils::find_first_above_threshold(
        tables[i].get_split_count_data(), n_tip, thresh);
      
      if (first_qualifying >= 0) {
        tables[i].find_splits_above_threshold_simd(thresh, qualifying_splits);
        for (int32 j = tables[i].X_left(k + 1);
             j != tables[i].X_right(k + 1) + 1;
             ++j) {
          ret(splits_found, tables[i].DECODE(j) - 1) = true;
        }
        ++splits_found;
        if (splits_found == ntip_3) {
          return ret;
        }
      }
    }
    
  } while (i++ != n_trees - thresh);
  
  return splits_found ? 
  ret(Range(0, splits_found - 1), _) :
    LogicalMatrix(0, n_tip);
}
