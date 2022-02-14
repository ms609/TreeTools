#include <Rcpp/Lightest>
#include <random>
#include "../inst/include/TreeTools.h"
using namespace Rcpp;

const intx MAX_TIP = 44, MAX_NODE = MAX_TIP + MAX_TIP - 1;

// [[Rcpp::export]]
IntegerVector num_to_parent(const NumericVector n, const IntegerVector nTip) {
  if (nTip[0] < 2) {
    throw std::range_error("nTip must be > 1");
  }
  const intx
    n_tip = nTip[0],
    root_node = n_tip + n_tip - 1,
    c_to_r = 1,
    prime = n_tip - 2
  ;
  uint64_t tree_id = n[0];
  for (intx i = 1; i < n.length(); i++) {
    tree_id *= INT_MAX;
    tree_id += n[i];
  }
  intx base;

  IntegerVector edge(n_tip + n_tip - 2);
  edge(0) = root_node;
  edge(1) = root_node;

  for (intx i = 2; i != n_tip; i++) {
    base = (i + i - 3);
    const intx
      i_prime = i + prime,
      i_prime_r = i_prime + c_to_r
    ;

    intx where = (tree_id % base) + 1;
    if (where >= i) {
      where += prime + 2 - i;
    }

    edge(i_prime) = edge(where);
    edge(i) = i_prime_r;
    edge(where) = i_prime_r;

    tree_id /= base;
  }

  return edge;
}

// [[Rcpp::export]]
IntegerVector random_parent(const IntegerVector nTip, const IntegerVector seed) {
  if (nTip[0] < 2) {
    throw std::range_error("nTip must be > 1");
  }
  const intx
    n_tip = nTip[0],
    root_node = n_tip + n_tip - 1,
    c_to_r = 1,
    prime = n_tip - 2
  ;
  intx base;
  
  std::mt19937 rng;
  rng.seed(seed[0]);

  IntegerVector edge(n_tip + n_tip - 2);
  edge(0) = root_node;
  edge(1) = root_node;

  for (intx i = 2; i != n_tip; i++) {
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

intx minimum (const intx x, const intx y) {
  return (x < y) ? x : y;
}

intx maximum (const intx x, const intx y) {
  return (x > y) ? x : y;
}

// Parent and child must be in postorder, with tree rooted on tip 1.
// [[Rcpp::export]]
IntegerVector edge_to_num(IntegerVector parent, IntegerVector child,
                   IntegerVector nTip) {
  if (parent.size() != child.size()) {
    throw std::length_error("Parent and child must be the same length");
  }
  const intx
    n_tip = nTip[0],
    n_internal = n_tip - 1,
    n_edge = parent.size(),
    all_node = n_internal + n_tip,
    r_to_c = 1
  ;
  if (n_tip < 4) {
    return IntegerVector(1); // A length one, zero initialized vector
  }
  if (n_edge != n_tip + n_tip - 2) {
    throw std::length_error("nEdge must == nTip + nTip - 2");
  }
  intx smallest_below[MAX_NODE],
       parent_of[MAX_NODE],
       prime_id[MAX_NODE],
       index[MAX_TIP]
  ;
  for (intx i = 0; i != all_node; i++) {
    smallest_below[i] = i;
    prime_id[i] = i;
  }

  for (intx i = 0; i != n_edge - 2; i += 2) {
    const intx
      this_node = parent[i] - r_to_c,
      left_child = child[i] - r_to_c,
      right_child = child[i + 1] - r_to_c
    ;
    smallest_below[this_node] = minimum(smallest_below[right_child],
                                        smallest_below[left_child]);
    prime_id[this_node] = maximum(smallest_below[left_child],
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
  uint64_t num = 0U;
  uint64_t multiplier = 1U;
  for (intx i = 3; i != n_tip; i++) {
    intx insertion_edge = index[i];
    if (insertion_edge < n_tip) {
      --insertion_edge;
    } else {
      insertion_edge += i - (n_tip + 3);
    }
    num += insertion_edge * multiplier;
    multiplier *= (i + i - 3);
  }

  if (num >= INT_MAX) {
    return IntegerVector{int(num / INT_MAX), int(num % INT_MAX)};
  } else {
    return IntegerVector{int(num)};
  }
}

// Parent and child must be in postorder, with tree rooted on tip 1.
// [[Rcpp::export]]
IntegerVector edge_to_mixed_base(IntegerVector parent, IntegerVector child,
                                 IntegerVector nTip) {
  if (parent.size() != child.size()) {
    throw std::length_error("Parent and child must be the same length");
  }
  const intx
    n_tip = nTip[0],
    n_internal = n_tip - 1,
    n_edge = parent.size(),
    all_node = n_internal + n_tip,
    r_to_c = 1
  ;
  if (n_tip < 4) {
    return IntegerVector(0);
  }
  if (n_edge != n_tip + n_tip - 2) {
    throw std::length_error("nEdge must == nTip + nTip - 2");
  }
  intx
    smallest_below[MAX_NODE],
    parent_of[MAX_NODE],
    prime_id[MAX_NODE],
    index[MAX_TIP]
  ;
  for (intx i = 0; i != all_node; i++) {
    smallest_below[i] = i;
    prime_id[i] = i;
  }
  
  for (intx i = 0; i != n_edge - 2; i += 2) {
    const intx
      this_node = parent[i] - r_to_c,
      left_child = child[i] - r_to_c,
      right_child = child[i + 1] - r_to_c
    ;
    smallest_below[this_node] = minimum(smallest_below[right_child],
                                        smallest_below[left_child]);
    prime_id[this_node] = maximum(smallest_below[left_child],
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
  IntegerVector ret(n_tip - 3);
  for (intx i = 3; i != n_tip; i++) {
    intx insertion_edge = index[i];
    if (insertion_edge < n_tip) {
      --insertion_edge;
    } else {
      insertion_edge += i - (n_tip + 3);
    }
    ret[n_tip - i - 1] = insertion_edge;
  }
  return ret;
}

// [[Rcpp::export]]
IntegerVector mixed_base_to_parent(
    const IntegerVector n,
    const IntegerVector nTip
  ) {
  if (nTip[0] < 2) {
    throw std::range_error("nTip must be > 1");
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
