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

// equivalent to consensus_tree
// [[Rcpp::export]]
List count_splits(const List trees) {
  int16
    L_i, R_i,
    L_j, R_j
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

  IntegerVector split_n (CT_MAX_LEAVES);
  NumericVector split_pi (CT_MAX_LEAVES);
  NumericVector split_ci (CT_MAX_LEAVES);
  LogicalMatrix split_members(CT_MAX_LEAVES, n_tip);

  int32
    i = 0,
    splits_found = 0
  ;
  do { // For i in 0 to n_trees-2
    Rcout << "Tree " << i << ": ==================\n";
    for (int32 ii = 0; ii < ntip_3; ii++) Rcout << tables[i].X(ii, 0) << ", "
    << tables[i].X(ii, 1) << "\n";
    if (tables[i].NOSWX(ntip_3)) {
      // All switches are already "on"; thus clades have been counted and 
      // there's nothing new to see
      continue;
    }
    tables[i].XRESET(); // Reset enumeration to zero
    for (int16 sp = 0; sp != ntip_3; sp++) {
      // Read next cluster into L_i and R_i
      tables[i].NCLUS(&L_i, &R_i);
      if (tables[i].GETSWX()) {
        Rcout << "Already counted. \n";
        // Already counted
        continue;
      }
      if (!L_i) {
        Rcout << "Invalid cluster. \n";
        continue; // Non-cluster
      }
      Rcout << L_i << "-" << R_i << " represents leaves " << 
        tables[i].DECODE(L_i) << "-" << tables[i].DECODE(R_i) << "\n";
      
      int32 in_split = R_i - L_i + 1;
      Rcout << in_split << " leaves\n";
      for (int32 j = L_i; j != R_i + 1; ++j) {
        Rcout << tables[i].DECODE(j) << ", ";
        split_members(splits_found, tables[i].DECODE(j) - 1) = true;
      }
      ++(split_n[splits_found]); // Count occurrence in tree i
      for (int32 j = i + 1; j != n_trees; j++) {
        L_j = tables[j].ENCODE(tables[i].DECODE(L_i));
        R_j = tables[j].ENCODE(tables[i].DECODE(R_i));
        if (tables[j].ISCLUST(&L_j, &R_j)) {
          Rcout << "Present in tree " << j << "\n";
          tables[j].SETSW(&L_j, &R_j);
          ++(split_n[splits_found]);
        }
      }
      Rcout << "Found " << split_n[splits_found] << " times. ";
      

      
      split_pi[splits_found] = TreeTools::split_phylo_info(
        int16(in_split), &n_tip_16, split_n[splits_found] / double(n_trees));
      
      // Rcout << "Phylo Information content: " << split_pi[splits_found] << "\n";
      
      split_ci[splits_found] = TreeTools::split_clust_info(
        int16(in_split), &n_tip_16, split_n[splits_found] / double(n_trees));
      
      // For completeness:
      // tables[i].SETSWX()
      // No need to call though: we'll never visit again.
      ++splits_found;
      Rcout << "\n That's " << splits_found << " in the bank!\n\n";
    }
  } while (++i != n_trees - 1);
  

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
List count_splits_bku(const List trees) {
  int16
    v = 0, w = 0,
    L, R, N, W,
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
            Rcout << "Split" << i << ", " << j << "'s been counted\n";
            // Split has already been counted; next!
          } else {
            if (N == R - L + 1) { // L..R is contiguous, and must be tested
              if (tables[i].CLUSTONL(&L, &R)) {
                tables[j].SETSWX(&j_pos);
                ASSERT(L > 0);
                Rcout << "Split " << i << ", "<<j<<" counting L(" << (L-1)
                  << ")\n";
                ++split_count[L - 1];
              } else if (tables[i].CLUSTONR(&L, &R)) {
                tables[j].SETSWX(&j_pos);
                ASSERT(R > 0);
                Rcout << "Split " << i << ", "<<j<<" counting R(" << (R-1)
                      << ")\n";
                ++split_count[R - 1];
              }
            }
          }
        }
        tables[j].NVERTEX_short(&v, &w);
      } while (v);
    }

    // Compute information content of each split
    for (int32 k = n_tip; k--; ) {
      int32 in_split = tables[i].X(k + 1, 1) - tables[i].X(k + 1, 0) + 1;
      if (in_split > 1 && in_split < n_tip - 1) {
        
        if (n_tip > CT_MAX_LEAVES) {
          Rcerr << "Too many leaves";
          return List::create();
        }
        
        split_n[splits_found] = split_count[k];
        
        int16 n_tip_16 = int16(n_tip);
        
        // Rcout << splits_found << ": Found tree " << i << "'s split " << k
        //       << " in " << split_count[k] << " trees.\n";
        // Rcout << "Count: " << split_count[k] <<"\n";
        
        for (int32 j = tables[i].X(k + 1, 0);
             j != tables[i].X(k + 1, 1) + 1;
             ++j) {
          // Rcout << ", " << tables[i].DECODE(j);
          split_members(splits_found, tables[i].DECODE(j) - 1) = true;
        }
        
        split_pi[splits_found] = TreeTools::split_phylo_info(
          int16(in_split), &n_tip_16, split_count[k] / double(n_trees));
        
        // Rcout << "Phylo Information content: " << split_pi[splits_found] << "\n";
        
        split_ci[splits_found] = TreeTools::split_clust_info(
          int16(in_split), &n_tip_16, split_count[k] / double(n_trees));
        
        // Rcout << "Clust Information content: " << split_ci[splits_found] << "\n";
        
        ++splits_found;
        // Rcout << "\n\n";
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
