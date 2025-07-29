#include <Rcpp/Lightest>
using namespace Rcpp;

#include "../inst/include/TreeTools/assert.h" /* for ASSERT */
#include "../inst/include/TreeTools/ClusterTable.h" /* for ClusterTable */
#include "../inst/include/TreeTools/ClusterTableSIMD.h" /* for SIMD extensions */

#include <algorithm> /* for fill */
#include <array> /* for array */
#include <vector> /* for vector */
#include <chrono> /* for benchmarking */

// Ultra-minimal test - just see if we can create ClusterTableSIMD
// [[Rcpp::export]]
bool test_basic_simd_construction(const List trees) {
  try {
    // Just test construction of first tree
    TreeTools::ClusterTable orig(trees[0]);
    TreeTools::ClusterTableSIMD simd_table(trees[0]);
    
    // Test that basic methods work
    int orig_n = orig.N();
    int simd_n = simd_table.N();
    
    return (orig_n == simd_n && orig_n > 0);
    
  } catch (...) {
    return false;
  }
}

// Test just the SIMD utilities without ClusterTable
// [[Rcpp::export]]
List test_pure_simd() {
  std::vector<int32> test_data = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  int32 threshold = 5;
  
  int32 simd_count = TreeTools::simd_utils::count_above_threshold(
    test_data.data(), test_data.size(), threshold);
  
  int32 manual_count = 0;
  for (int32 val : test_data) {
    if (val >= threshold) ++manual_count;
  }
  
  return List::create(
    Named("simd_count") = simd_count,
    Named("manual_count") = manual_count,
    Named("match") = (simd_count == manual_count)
  );
}

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

// Simple benchmark function with better error handling
// [[Rcpp::export]]
List benchmark_consensus_simple(const List trees, const NumericVector p, 
                                int n_iterations = 10) {
  
  LogicalMatrix result_original, result_simd;
  double time_original, time_simd;
  bool results_match = false;
  std::string error_msg = "";
  
  try {
    // Test original function
    auto start_original = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_iterations; ++i) {
      result_original = consensus_tree(trees, p);
    }
    auto end_original = std::chrono::high_resolution_clock::now();
    time_original = std::chrono::duration<double, std::milli>(
      end_original - start_original).count();
    
    // Test SIMD function
    auto start_simd = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_iterations; ++i) {
      result_simd = consensus_tree_simd_minimal(trees, p);
    }
    auto end_simd = std::chrono::high_resolution_clock::now();
    time_simd = std::chrono::duration<double, std::milli>(
      end_simd - start_simd).count();
    
    // Check dimensions first
    if (result_original.nrow() != result_simd.nrow() || 
        result_original.ncol() != result_simd.ncol()) {
      results_match = false;
      error_msg = "Dimension mismatch: Original(" + 
        std::to_string(result_original.nrow()) + "x" + 
        std::to_string(result_original.ncol()) + ") vs SIMD(" +
        std::to_string(result_simd.nrow()) + "x" + 
        std::to_string(result_simd.ncol()) + ")";
    } else {
      // Compare element by element
      results_match = true;
      for (int i = 0; i < result_original.nrow() && results_match; ++i) {
        for (int j = 0; j < result_original.ncol() && results_match; ++j) {
          if (result_original(i, j) != result_simd(i, j)) {
            results_match = false;
            error_msg = "Element mismatch at (" + std::to_string(i) + "," + 
              std::to_string(j) + "): " + 
              std::to_string(result_original(i, j)) + " vs " + 
              std::to_string(result_simd(i, j));
          }
        }
      }
    }
    
  } catch (const std::exception& e) {
    error_msg = std::string("Exception: ") + e.what();
    time_original = -1;
    time_simd = -1;
  }
  
  return List::create(
    Named("time_original_ms") = time_original,
    Named("time_simd_ms") = time_simd,
    Named("speedup") = time_simd > 0 ? time_original / time_simd : -1,
    Named("results_match") = results_match,
    Named("error_message") = error_msg,
    Named("original_dims") = IntegerVector::create(result_original.nrow(), result_original.ncol()),
    Named("simd_dims") = IntegerVector::create(result_simd.nrow(), result_simd.ncol()),
    Named("simd_features") = List::create(
      Named("avx2") = TreeTools::simd_utils::has_avx2_support(),
      Named("sse2") = TreeTools::simd_utils::has_sse2_support()
    )
  );
}

// Just test if we can run consensus_tree with regular ClusterTable
// [[Rcpp::export]]
bool test_original_consensus(const List trees, const NumericVector p) {
  try {
    LogicalMatrix result = consensus_tree(trees, p);
    return true;
  } catch (...) {
    return false;
  }
}