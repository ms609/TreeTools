#include <Rcpp/Lightest>
using namespace Rcpp;

#include "../inst/include/TreeTools/assert.h" /* for ASSERT */
#include "../inst/include/TreeTools/ClusterTable.h" /* for ClusterTable */

#include <algorithm> /* for fill */
#include <array> /* for array */
#include <vector> /* for vector */

// trees is a list of objects of class phylo, all with the same tip labels
// (try RenumberTips(trees, trees[[1]]))
// Per #168, unexpected behaviour if root position differs in non-preorder trees
// Further investigation could be beneficial; for now, suggest applying
// the function to preorder trees only.
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

  int32
    i = 0,
    splits_found = 0
  ;
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
        // Rcout << splits_found << ": Found tree " << i << "'s split " << k
        //       << " in " << split_count[k] << " trees.\n";
        for (int32 j = tables[i].X(k + 1, 0);
             j != tables[i].X(k + 1, 1) + 1;
             ++j) {
          // Rcout << ", " << tables[i].DECODE(j);
          ret(splits_found, tables[i].DECODE(j) - 1) = true;
        }
        // Rcout << "\n\n";
        ++splits_found;
        // If we have a perfectly resolved tree, break.
        if (splits_found == ntip_3) {
          return ret;
        }
      }
    }
  } while (i++ != n_trees - thresh); // All clades in p% consensus must occur in first q% of trees.

  return splits_found ? 
    ret(Range(0, splits_found - 1), _) :
    LogicalMatrix(0, n_tip);
}
