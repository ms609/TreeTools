#include <Rcpp/Lighter> /* for is_na */
#include <algorithm>    /* for std::reverse, std::min, std::max */
#include <random>
#include <vector>
#include "../inst/include/TreeTools.h"
#include "../inst/include/TreeTools/tree_number.h"
using namespace Rcpp;

constexpr intx MB_MAX_TIP = 32768;
constexpr intx MB_MAX_NODE = MB_MAX_TIP + MB_MAX_TIP - 1;

// Build a tree_num_t from an INT_MAX-packed IntegerVector.
// Encoding: tree_id = n[0] * INT_MAX^(k-1) + ... + n[k-1]  (most-significant first)
inline TreeTools::tree_num_t packed_to_tree_num(const IntegerVector& n) {
  TreeTools::tree_num_t result(static_cast<uint64_t>(n[0]));
  for (int i = 1; i < n.length(); ++i) {
    result.mul_small(static_cast<uint64_t>(INT_MAX));
    result.add_small(static_cast<uint64_t>(n[i]));
  }
  return result;
}

// Convert a tree_num_t to an INT_MAX-packed IntegerVector.
// Digits are stored most-significant first, matching packed_to_tree_num.
inline IntegerVector tree_num_to_packed(TreeTools::tree_num_t num) {
  const TreeTools::tree_num_t zero;
  if (num == zero) {
    return IntegerVector{0};
  }
  std::vector<int> digits;
  while (!(num == zero)) {
    digits.push_back(
      static_cast<int>(num.divmod_small(static_cast<uint64_t>(INT_MAX))));
  }
  std::reverse(digits.begin(), digits.end());
  return IntegerVector(digits.begin(), digits.end());
}

// [[Rcpp::export]]
IntegerVector num_to_parent(const IntegerVector n, const IntegerVector nTip) {
  if (Rcpp::is_true(Rcpp::any(Rcpp::is_na(n)))) {
    Rcpp::stop("`n` may not contain NA values");
  }
  if (Rcpp::is_true(Rcpp::any(n < 0))) {
    Rcpp::stop("`n` may not be negative");
  }
  if (nTip[0] < 2) {
    Rcpp::stop("`nTip` must be > 1");
  }
  if (nTip.length() > 1) {
    Rcpp::warning("`nTip` should be a single integer");
  }
  const intx n_tip = nTip[0];
  const intx n_edge = n_tip + n_tip - 2;

  const TreeTools::tree_num_t tree_id = packed_to_tree_num(n);

  std::vector<intx> buf(n_edge);
  TreeTools::tree_number_to_parent(tree_id, n_tip, buf.data());

  IntegerVector edge(n_edge);
  for (intx i = 0; i < n_edge; ++i) {
    edge[i] = static_cast<int>(buf[i]);
  }
  return edge;
}

// Checking that nTip > 2 is caller's responsibility.
// [[Rcpp::export]]
IntegerVector random_parent(const IntegerVector nTip, const IntegerVector seed) {

  const intx
    n_tip = nTip[0],
    root_node = n_tip + n_tip - 1,
    c_to_r = 1,
    prime = n_tip - 2
  ;
  intx base;

  std::mt19937 rng(seed[0]);

  IntegerVector edge(n_tip + n_tip - 2);
  edge(0) = root_node;
  edge(1) = root_node;
  
  const intx i2 = 2;
  const intx where2 = 1;
  edge(i2 + prime) = edge(where2);
  edge(i2) = i2 + prime + c_to_r;
  edge(where2) = i2 + prime + c_to_r;

  for (intx i = 3; i != n_tip; i++) {
    base = (i + i - 3);
    const intx
      i_prime = i + prime,
      i_prime_r = i_prime + c_to_r
    ;
    
    std::uniform_int_distribution<std::mt19937::result_type> place(1, base);
    intx where = place(rng);
    if (where >= i) {
      where += prime + 2 - i;
    }

    edge(i_prime) = edge(where);
    edge(i) = i_prime_r;
    edge(where) = i_prime_r;
  }

  return edge;
}

// Parent and child must be in postorder, with tree rooted on tip 1.
// [[Rcpp::export]]
IntegerVector edge_to_num(
    const IntegerVector& parent,
    const IntegerVector& child,
    const IntegerVector& nTip) {
  if (parent.size() != child.size()) {
    Rcpp::stop("Parent and child must be the same length");
  }
  const intx n_tip = nTip[0];
  const intx n_edge = parent.size();

  if (n_tip < 4) {
    return IntegerVector(1); // A length one, zero initialized vector
  }
  if (n_edge != n_tip + n_tip - 2) {
    Rcpp::stop("nEdge must == nTip + nTip - 2");
  }
  if (n_tip > TreeTools::TREE_NUM_MAX_TIP) {
    Rcpp::stop("Too many leaves for tree number representation");
  }

  std::vector<intx> par_arr(n_edge), chi_arr(n_edge);
  for (intx i = 0; i < n_edge; ++i) {
    par_arr[i] = parent[i];
    chi_arr[i] = child[i];
  }

  return tree_num_to_packed(
    TreeTools::edges_to_tree_number(par_arr.data(), chi_arr.data(), n_tip));
}

inline void calc_edge_to_mixed_base(
    const IntegerVector& parent,
    const IntegerVector& child,
    const intx n_tip,
    const intx all_node,
    const intx n_edge,
    IntegerVector& ret) {
  
  constexpr intx r_to_c = 1;
  
  intx smallest_below[MB_MAX_NODE];
  intx parent_of[MB_MAX_NODE];
  intx prime_id[MB_MAX_NODE];
  intx index[MB_MAX_TIP];
  
  for (intx i = 0; i != all_node; i++) {
    smallest_below[i] = i;
    prime_id[i] = i;
  }
  
  for (intx i = 0; i < n_edge - 2; i += 2) {
    const intx this_node = parent[i] - r_to_c;
    const intx left_child = child[i] - r_to_c;
    const intx right_child = child[i + 1] - r_to_c;
    
    smallest_below[this_node] = std::min(smallest_below[right_child],
                                         smallest_below[left_child]);
    prime_id[this_node] = std::max(smallest_below[left_child],
                                   smallest_below[right_child]);
    parent_of[left_child] = parent_of[right_child] = this_node;
    
    for (intx at = smallest_below[this_node];
         at != this_node;
         at = parent_of[at]) {
      const intx prime_candidate = prime_id[at];
      if (prime_candidate < prime_id[this_node]) {
        index[prime_id[this_node]] = prime_id[at] + (at < n_tip ? 0 : n_tip);
      }
    }
  }
  for (intx i = 3; i < n_tip; ++i) {
    intx insertion_edge = index[i];
    if (insertion_edge < n_tip) {
      --insertion_edge;
    } else {
      insertion_edge += i - (n_tip + 3);
    }
    ret[n_tip - i - 1] = insertion_edge;
  }
}

// Parent and child must be in postorder, with tree rooted on tip 1.
// [[Rcpp::export]]
IntegerVector edge_to_mixed_base(
    const IntegerVector parent,
    const IntegerVector child,
    const IntegerVector nTip) {
  if (parent.size() != child.size()) {
    Rcpp::stop("Parent and child must be the same length");
  }
  if (nTip.length() > 1) {
    Rcpp::warning("`nTip` should be a single integer");
  }
  const intx n_tip = nTip[0];
  const intx n_internal = n_tip - 1;
  const intx n_edge = parent.size();
  const intx all_node = n_internal + n_tip;
  
  if (n_tip < 4) {
    return IntegerVector(0);
  }
  if (n_edge != n_tip + n_tip - 2) {
    Rcpp::stop("nEdge must == nTip + nTip - 2");
  }
  if (all_node > MB_MAX_NODE) {
    Rcpp::stop("Too many nodes for mixed base representation");
  }
  if (n_tip >= MB_MAX_TIP) {
    Rcpp::stop("Too many leaves for mixed base representation");
  }
  IntegerVector ret(n_tip - 3);
  calc_edge_to_mixed_base(parent, child, n_tip, all_node, n_edge, ret);
  return ret;
}

// [[Rcpp::export]]
IntegerVector mixed_base_to_parent(
    const IntegerVector n,
    const IntegerVector nTip
  ) {
  if (Rcpp::is_true(Rcpp::any(Rcpp::is_na(n)))) {
    Rcpp::stop("`n` may not contain NA values");
  }
  if (Rcpp::is_true(Rcpp::any(n < 0))) {
    Rcpp::stop("`n` may not be negative");
  }
  if (nTip[0] < 2) {
    Rcpp::stop("`nTip` must be > 1");
  }
  if (nTip.length() > 1) {
    Rcpp::warning("`nTip` should be a single integer");
  }
  const intx
    n_tip = nTip[0],
    root_node = n_tip + n_tip - 1,
    c_to_r = 1,
    prime = n_tip - 2
  ;
  
  IntegerVector edge(n_tip + n_tip - 2);
  edge(0) = root_node;
  edge(1) = root_node;
  
  for (intx i = 2; i != n_tip; i++) {
    const intx
      i_prime = i + prime,
      i_prime_r = i_prime + c_to_r
    ;
    
    intx where = i == 2 ? 1 : n[n_tip - i - 1] + 1;
    if (where >= i) {
      where += prime + 2 - i;
    }
    
    edge(i_prime) = edge(where);
    edge(i) = i_prime_r;
    edge(where) = i_prime_r;
  }
  
  return edge;
}
