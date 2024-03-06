#include <Rcpp/Lightest>
using namespace Rcpp;

#include "../inst/include/TreeTools/assert.h" /* for ASSERT */
#include "../inst/include/TreeTools/ClusterTable.h" /* for ClusterTable */
#include "../inst/include/TreeTools/information.h" /* for split_xxx_info */

#include <algorithm> /* for fill */
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

    // TODO Check: should this be n_tip - 1; or stop at k = 1 instead of 0?
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

// TODO document parameters and return value
// [[Rcpp::export]]
List count_splits(const List trees) {
  int16
    v = 0, w = 0,
    L, R, N, W,
    L_i, R_i, count,
    L_j, R_j, N_j, W_j
  ;
  const int32 n_trees = trees.length();

  std::vector<TreeTools::ClusterTable> tables;
  tables.reserve(n_trees);
  for (int32 i = n_trees; i--; ) {
    tables.emplace_back(TreeTools::ClusterTable(Rcpp::List(trees(i))));
  }

  const int32
    n_tip = tables[0].N(),
    ntip_3 = n_tip - 3
  ;
  const int16 n_tip_16 = int16(n_tip);

  std::array<int32, CT_STACK_SIZE * CT_MAX_LEAVES> S;
  std::array<int32, CT_MAX_LEAVES> split_count;
  
  IntegerVector split_n (ntip_3 * n_trees);
  NumericVector split_pi (ntip_3 * n_trees);
  NumericVector split_ci (ntip_3 * n_trees);
  LogicalMatrix split_members(ntip_3 * n_trees, n_tip);

  int32
    i = 0,
    splits_found = 0
  ;
  do {
    std::fill(split_count.begin(), split_count.begin() + n_tip, 1);

    for (int32 j = 0; j != n_trees; j++) {
      if (i == j) continue;

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
          
          if (N == R - L + 1) { // L..R is contiguous, and must be tested
            if (tables[i].CLUSTONL(&L, &R)) {
              ASSERT(L > 0);
              ++split_count[L - 1];
            } else if (tables[i].CLUSTONR(&L, &R)) {
              ASSERT(R > 0);
              ++split_count[R - 1];
            }
          
          }
        }
        tables[j].NVERTEX_short(&v, &w);
      } while (v);
    }
    
    // TODO Check: should this be n_tip - 1; or stop at k = 1 instead of 0?
    for (int32 k = 1; k <= ntip_3; k++) {
      L_i = tables[i].X(k + 1, 0);
      R_i = tables[i].X(k + 1, 1);
      count = split_count[k];
      Rcout << L_i << "-"<<R_i<<"...\n";
      int32 in_split = R_i - L_i + 1;
      if (count &&
          in_split > 1 &&
          in_split < n_tip - 1
      ) {
        Rcout << splits_found << ": Found tree " << i << "'s split " << k
              << " in " << count << " trees.\n";
        for (int32 leaf = L_i; leaf != R_i + 1; ++leaf) {
          Rcout << ", " << tables[i].DECODE(leaf);
          split_members(splits_found, tables[i].DECODE(leaf) - 1) = true;
        }
        
        split_n[splits_found] = count;
        
        split_pi[splits_found] = TreeTools::split_phylo_info(
          int16(in_split), &n_tip_16, split_n[splits_found] / double(n_trees));
        
        split_ci[splits_found] = TreeTools::split_clust_info(
          int16(in_split), &n_tip_16, split_n[splits_found] / double(n_trees));
        
        Rcout << "\n\n";
        ++splits_found;
        
      }
    }
  } while (++i != n_trees);

  return splits_found ? 
    List::create(
      Named("splits") = split_members(Range(0, splits_found - 1), _),
      _["count"] = split_n[Range(0, splits_found - 1)],
      _["pic"] = split_pi[Range(0, splits_found - 1)],
      _["cic"] = split_ci[Range(0, splits_found - 1)]
  ) : List::create(
      Named("splits") = LogicalMatrix(0, n_tip),
      _["count"] = IntegerVector::create(0),
      _["pic"] = NumericVector::create(0),
      _["cic"] = NumericVector::create(0)
  );
}

// equivalent to consensus_tree
// [[Rcpp::export]]
List count_splits_wrongly(const List trees) {
  int16
    v = 0, w = 0,
    L, R, N, W,
    L_i, R_i,
    L_j, R_j, N_j, W_j,
    count
  ;
  const int32 n_trees = trees.length();
  
  std::vector<TreeTools::ClusterTable> tables;
  tables.reserve(n_trees);
  for (int32 i = n_trees; i--; ) {
    tables.emplace_back(TreeTools::ClusterTable(Rcpp::List(trees(i))));
  }

  const int32
    n_tip = tables[0].N(),
    ntip_3 = n_tip - 3
  ;
  const int16 n_tip_16 = int16(n_tip);

  std::array<int32, CT_STACK_SIZE * CT_MAX_LEAVES> S;
  std::array<int32, CT_MAX_LEAVES> split_count;
  
  IntegerVector split_n (CT_MAX_LEAVES);
  NumericVector split_pi (CT_MAX_LEAVES);
  NumericVector split_ci (CT_MAX_LEAVES);
  LogicalMatrix split_members(CT_MAX_LEAVES, n_tip);

  int32
    i = 0,
    splits_found = 0
  ;
  do {
    if (tables[i].NOSWX(ntip_3)) {
      continue;
    }
    Rcout << "\nTree " << i << " ===========\n";

    // UNTESTED ASSERTION: we can change to n_tip - 1, as the final entry is
    // all tips, thus is uninteresting.
    for (int16 it = 1; it != n_tip - 1; it++) {
      split_count[it - 1] = 1 - tables[i].GETSWX(&it);
    }
    for (int16 j = 0; j != n_tip - 2; ++j) {
      Rcout << split_count[j];
    }
    Rcout << "\n\n";

    for (int32 j = i + 1; j != n_trees; j++) {

      //tables[i].CLEAR();
      Rcout << "\n * Reading table " << j << "; " << tables[j].NOSWX() << 
        " switches on.\n";

      tables[j].TRESET();
      tables[j].READT(&v, &w);

      int16 j_pos = 0;
      int16 j_row = 0;
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
          // Rcout << " j Row " << j_pos << ": " << tables[j].X(j_pos, 0)
          //   << "-" << tables[j].X(j_pos, 1) << "\n";
          Rcout << " Examining I's " << L << "-" << R<<"\n";
          j_row = (tables[j].CLUSTONL(&L_j, &R_j) ? L_j : R_j);
          Rcout << " Row " << j_row << " on J: " << tables[j].X(j_row, 0)
            << "-" << tables[j].X(j_row, 1) << "\n";
          if (tables[j].GETSWX(&j_row)) {
            // Split has already been counted; next!
            Rcout << "Already counted at row " << j_row << ".\n";
          } else {
            if (N == R - L + 1) { // L..R is contiguous, and must be tested
              if (tables[i].CLUSTONL(&L, &R)) {
                Rcout << "Counting " << L << "-" << R << " on tree " << j
                << " at pos " << j_pos << "L.\n";
                tables[j].SETSWX(&j_row);
                
                ASSERT(L > 0);
                ++split_count[L - 1];
              } else if (tables[i].CLUSTONR(&L, &R)) {
                Rcout << "Counting tree " << j << " at pos " << j_pos << "R.\n";
                tables[j].SETSWX(&j_row);
                ASSERT(R > 0);
                ++split_count[R - 1];
              }
            }
          }
        }
        tables[j].NVERTEX_short(&v, &w);
      } while (v);
    }

    for (int32 k = n_tip - 2; k--; ) {
      L_i = tables[i].X(k + 1, 0);
      R_i = tables[i].X(k + 1, 1);
      count = split_count[k];
      Rcout << L_i << "-"<<R_i<<"...\n";
      int32 in_split = R_i - L_i + 1;
      if (//count &&
          in_split > 1 &&
          in_split < n_tip - 1
        ) {
        Rcout << splits_found << ": Found tree " << i << "'s split " << k
              << " in " << count << " trees.\n";
        for (int32 leaf = L_i; leaf != R_i + 1; ++leaf) {
          Rcout << ", " << tables[i].DECODE(leaf);
          if (count)//////////////////
          split_members(splits_found, tables[i].DECODE(leaf) - 1) = true;
        }
        if (!count) continue; ////////////////
        
        split_n[splits_found] = count;
          
        split_pi[splits_found] = TreeTools::split_phylo_info(
          int16(in_split), &n_tip_16, split_n[splits_found] / double(n_trees));
        
        split_ci[splits_found] = TreeTools::split_clust_info(
          int16(in_split), &n_tip_16, split_n[splits_found] / double(n_trees));
        
        Rcout << "\n\n";
        ++splits_found;
        
      }
    }
  } while (++i != n_trees);

  Rcout << "Found " << splits_found << " splits in total. Returning. \n";
  
  return splits_found ? 
    List::create(
      Named("splits") = split_members(Range(0, splits_found - 1), _),
      _["count"] = split_n[Range(0, splits_found - 1)],
      _["pic"] = split_pi[Range(0, splits_found - 1)],
      _["cic"] = split_ci[Range(0, splits_found - 1)]
  ) : List::create(
      Named("splits") = LogicalMatrix(0, n_tip),
      _["count"] = IntegerVector::create(0),
      _["pic"] = NumericVector::create(0),
      _["cic"] = NumericVector::create(0)
  );
}
