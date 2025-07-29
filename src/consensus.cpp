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

// ULTRA-MINIMAL SIMD VERSION - Just test if ClusterTableSIMD works in the loop
// [[Rcpp::export]]
LogicalMatrix consensus_tree_simd_minimal(const List trees, const NumericVector p) {
  // First, let's just test if we can create ClusterTableSIMD objects and run original logic
  
  const int32
  n_trees = trees.length(),
    frac_thresh = int32(n_trees * p[0]) + 1,
    thresh = frac_thresh > n_trees ? n_trees : frac_thresh
  ;
  
  // Try to use SIMD-enabled ClusterTable but with original algorithm
  std::vector<TreeTools::ClusterTableSIMD> tables;
  tables.reserve(n_trees);
  
  // Construct each table carefully
  for (int32 i = 0; i < n_trees; ++i) {
    const Rcpp::List& tree_i = Rcpp::List(trees(i));
    tables.emplace_back(TreeTools::ClusterTableSIMD(tree_i));
    
    // Verify construction worked
    if (!tables[i].is_properly_initialized()) {
      Rcpp::stop("ClusterTableSIMD construction failed for tree " + std::to_string(i));
    }
  }
  
  const int32
  n_tip = tables[0].N(),
    ntip_3 = n_tip - 3
  ;
  
  // Use original algorithm but with SIMD tables
  std::array<int32, CT_STACK_SIZE * CT_MAX_LEAVES> S;
  std::array<int32, CT_MAX_LEAVES> split_count;  // Keep original split_count for now
  
  LogicalMatrix ret(ntip_3, n_tip);
  
  int32 i = 0;
  int32 splits_found = 0;
  
  int16 v = 0, w = 0, L, R, N, W, L_j, R_j, N_j, W_j;
  
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
                ++split_count[L - 1];  // Keep original logic
              } else if (tables[i].CLUSTONR(&L, &R)) {
                tables[j].SETSWX(&j_pos);
                ASSERT(R > 0);
                ++split_count[R - 1];  // Keep original logic
              }
            }
          }
        }
        tables[j].NVERTEX_short(&v, &w);
      } while (v);
    }
    
    // Original threshold checking
    for (int32 k = n_tip; k--; ) {
      if (split_count[k] >= thresh) {
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

// Version 2: Add actual SIMD split counting
// [[Rcpp::export]]
LogicalMatrix consensus_tree_simd_v2(const List trees, const NumericVector p) {
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
    
    // NOW USE SIMD split counting instead of std::array
    tables[i].init_split_count_simd_aligned(n_tip);
    
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
                // USE SIMD increment
                tables[i].increment_split_count_simd(L - 1);
              } else if (tables[i].CLUSTONR(&L, &R)) {
                tables[j].SETSWX(&j_pos);
                ASSERT(R > 0);
                // USE SIMD increment
                tables[i].increment_split_count_simd(R - 1);
              }
            }
          }
        }
        tables[j].NVERTEX_short(&v, &w);
      } while (v);
    }
    
    // USE SIMD threshold checking
    for (int32 k = n_tip; k--; ) {
      if (tables[i].get_split_count_simd(k) >= thresh) {
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
// Version 3: Vectorized threshold checking
// [[Rcpp::export]]
LogicalMatrix consensus_tree_simd_v3(const List trees, const NumericVector p) {
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
    
    tables[i].init_split_count_simd_aligned(n_tip);
    
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
                tables[i].increment_split_count_simd_aligned(L - 1);
              } else if (tables[i].CLUSTONR(&L, &R)) {
                tables[j].SETSWX(&j_pos);
                ASSERT(R > 0);
                tables[i].increment_split_count_simd_aligned(R - 1);
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

// Updated benchmark to test v2
// [[Rcpp::export]]
List benchmark_consensus_v2(const List trees, const NumericVector p, 
                            int n_iterations = 10) {
  
  LogicalMatrix result_original, result_simd_v1, result_simd_v2;
  double time_original, time_simd_v1, time_simd_v2;
  
  try {
    // Original
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_iterations; ++i) {
      result_original = consensus_tree(trees, p);
    }
    auto end = std::chrono::high_resolution_clock::now();
    time_original = std::chrono::duration<double, std::milli>(end - start).count();
    
    // SIMD v1 (minimal changes)
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_iterations; ++i) {
      result_simd_v1 = consensus_tree_simd_v3(trees, p);
    }
    end = std::chrono::high_resolution_clock::now();
    time_simd_v1 = std::chrono::duration<double, std::milli>(end - start).count();
    
    // SIMD v2 (real SIMD split counting)
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_iterations; ++i) {
      result_simd_v2 = consensus_tree_simd_v2(trees, p);
    }
    end = std::chrono::high_resolution_clock::now();
    time_simd_v2 = std::chrono::duration<double, std::milli>(end - start).count();
    
    // Check correctness
    bool v1_matches = (result_original.nrow() == result_simd_v1.nrow() && 
                       result_original.ncol() == result_simd_v1.ncol());
    bool v2_matches = (result_original.nrow() == result_simd_v2.nrow() && 
                       result_original.ncol() == result_simd_v2.ncol());
    
    if (v1_matches) {
      for (int i = 0; i < result_original.nrow() && v1_matches; ++i) {
        for (int j = 0; j < result_original.ncol() && v1_matches; ++j) {
          if (result_original(i, j) != result_simd_v1(i, j)) v1_matches = false;
        }
      }
    }
    
    if (v2_matches) {
      for (int i = 0; i < result_original.nrow() && v2_matches; ++i) {
        for (int j = 0; j < result_original.ncol() && v2_matches; ++j) {
          if (result_original(i, j) != result_simd_v2(i, j)) v2_matches = false;
        }
      }
    }
    
    return List::create(
      Named("time_original_ms") = time_original,
      Named("time_simd_v3_ms") = time_simd_v1,
      Named("time_simd_v2_ms") = time_simd_v2,
      Named("speedup_v1") = time_original / time_simd_v1,
      Named("speedup_v2") = time_original / time_simd_v2,
      Named("v1_matches") = v1_matches,
      Named("v2_matches") = v2_matches,
      Named("simd_features") = List::create(
        Named("avx2") = TreeTools::simd_utils::has_avx2_support(),
        Named("sse2") = TreeTools::simd_utils::has_sse2_support()
      )
    );
    
  } catch (const std::exception& e) {
    return List::create(
      Named("error") = std::string(e.what())
    );
  }
}
