#include <Rcpp/Lightest>
using namespace Rcpp;

#include "../inst/include/TreeTools/assert.h" /* for ASSERT */
#include "../inst/include/TreeTools/ClusterTable.h" /* for ClusterTable */

#include <algorithm> /* for fill */
#include <array> /* for array */

using TreeTools::ct_stack_threshold;
using TreeTools::ct_max_leaves_heap;

struct StackEntry { int32 L, R, N, W; };

// Helper template function to perform consensus computation
// Uses StackContainer for the S array (either std::array or std::vector)
template<typename StackContainer>
RawMatrix calc_consensus_tree(
  const List& trees,
  const NumericVector& p,
  StackContainer& S
) {
  int32 v = 0;
  int32 w = 0;
  int32 L, R, N, W;
  
  const int32 n_trees = trees.length();
  const int32 frac_thresh = int32(n_trees * p[0]) + 1;
  const int32 thresh = frac_thresh > n_trees ? n_trees : frac_thresh;
  
  std::vector<TreeTools::ClusterTable> tables;
  tables.reserve(n_trees);
  for (int32 i = 0; i < n_trees; ++i) {
    tables.emplace_back(TreeTools::ClusterTable(Rcpp::List(trees(i))));
  }
  
  const int32 n_tip   = tables[0].N();
  const int32 ntip_3  = n_tip - 3;
  const int32 nbin    = (n_tip + 7) / 8;  // bytes per row in packed output
  
  int32* split_count;
  std::array<int32, ct_stack_threshold> split_stack;
  std::vector<int32> split_heap;
  if (n_tip <= ct_stack_threshold) {
    split_count = split_stack.data();
  } else {
    split_heap.resize(n_tip);
    split_count = split_heap.data();
  }

  StackEntry *const S_start = S.data();
  
  // Packed output: each row has nbin bytes
  RawMatrix ret(ntip_3, nbin);

  int32 i = 0;
  int32 splits_found = 0;
  
  do {
    if (tables[i].NOSWX(ntip_3)) {
      continue;
    }
    
    std::fill(split_count, split_count + n_tip, 1);

    for (int32 j = i + 1; j < n_trees; ++j) {
      ASSERT(tables[i].N() == tables[j].N());

      tables[i].CLEAR();

      tables[j].TRESET();
      tables[j].READT(&v, &w);
      
      int32 j_pos = 0;
      StackEntry* S_top = S_start; // Empty the stack S
      
      do {
        if (CT_IS_LEAF(v)) {
          const auto enc_v = tables[i].ENCODE(v);
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
          
          if (!tables[j].GETSWX(&j_pos)) {
            if (N == R - L + 1) { // L..R is contiguous, and must be tested
              if (tables[i].CLUSTONL(L, R)) {
                tables[j].SETSWX(j_pos);
                ASSERT(L > 0);
                ++split_count[L - 1];
              } else if (tables[i].CLUSTONR(L, R)) {
                tables[j].SETSWX(j_pos);
                ASSERT(R > 0);
                ++split_count[R - 1];
              }
            }
          }
        }
        tables[j].NVERTEX_short(&v, &w);
      } while (v);
    }
    
    for (int32 k = 0; k < n_tip; ++k) {
      if (split_count[k] >= thresh) {
        const int32 start = tables[i].X_left(k + 1);
        const int32 end   = tables[i].X_right(k + 1);
        
        for (int32 j = start; j <= end; ++j) {
          const int32 leaf_idx = tables[i].DECODE(j) - 1; // 0-based
          const int32 byte_idx = leaf_idx >> 3;           // column index
          const int32 bit_idx  = leaf_idx & 7;            // bit within byte
          
          // pointer to the first row of this column
          Rbyte* col_ptr = &ret(0, byte_idx);
          col_ptr[splits_found] |= (Rbyte(1) << bit_idx); // set bit in row
        }
        
        ++splits_found;
        // If we have a perfectly resolved tree, exit early.
        if (splits_found == ntip_3) {
          return ret;
        }
      }
    }
  } while (i++ != n_trees - thresh); // All clades in p% consensus must occur in first q% of trees.
  
  return (splits_found == 0) ? RawMatrix(0, nbin) : 
         (splits_found < ntip_3) ? ret(Range(0, splits_found - 1), _) : ret;
}

// trees is a list of objects of class phylo, all with the same tip labels
// (try RenumberTips(trees, trees[[1]]))
// Per #168, unexpected behaviour if root position differs in non-preorder trees
// Further investigation could be beneficial; for now, suggest applying
// the function to preorder trees only.
// [[Rcpp::export]]
RawMatrix consensus_tree(const List trees, const NumericVector p) {
  // First, peek at the tree size to determine allocation strategy
  // We'll create a temporary ClusterTable just to check the size
  try {
    TreeTools::ClusterTable temp_table(Rcpp::List(trees(0)));
    const int32 n_tip = temp_table.N();
    
    if (n_tip <= ct_stack_threshold) {
      // Small tree: use stack-allocated array
      std::array<StackEntry, ct_stack_threshold> S;
      return calc_consensus_tree(trees, p, S);
    } else {
      // Large tree: use heap-allocated vector
      std::vector<StackEntry> S(n_tip);
      return calc_consensus_tree(trees, p, S);
    }
  } catch(const std::exception& e) {
    Rcpp::stop(e.what());
  }
  
  ASSERT(false && "Unreachable code in consensus_tree");
  return RawMatrix(0, 0);
}
