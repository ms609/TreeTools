#include <cstdint>
#include <Rcpp.h>
using namespace Rcpp;

const unsigned int MAX_TIP = 44, MAX_NODE = MAX_TIP + MAX_TIP - 1;
const uint64_t MAX_INT = 1000000000U;

// [[Rcpp::export]]
IntegerVector num_to_parent(NumericVector n, IntegerVector nTip) {
  if (nTip[0] < 2) {
    throw std::range_error("nTip must be > 1");
  }
  const unsigned int n_tip = nTip[0],
                     root_node = n_tip + n_tip - 1,
                     c_to_r = 1,
                     prime = n_tip - 2;
  uint64_t tree_id = n[0];
  for (int i = 1; i < n.length(); i++) {
    tree_id *= MAX_INT;
    tree_id += n[i];
  }
  unsigned int base;

  IntegerVector edge(n_tip + n_tip - 2);
  edge(0) = root_node;
  edge(1) = root_node;

  for (unsigned int i = 2; i < n_tip; i++) {
    base = (i + i - 3);
    const int i_prime = i + prime,
      i_prime_r = i_prime + c_to_r;

    uint64_t where = (tree_id % base) + 1;
    if (where >= i) {
      where += prime + 2 - i;
    }

    edge(i_prime) = edge(where);
    edge(i) = i_prime_r;
    edge(where) = i_prime_r;

    tree_id /= base;
  }

  return(edge);
}

unsigned int minimum (unsigned int x, unsigned int y) {
  return((x < y) ? x : y);
}

unsigned int maximum (unsigned int x, unsigned int y) {
  return((x > y) ? x : y);
}

// Parent and child must be in postorder, with tree rooted on tip 1.
// [[Rcpp::export]]
NumericVector edge_to_num(IntegerVector parent, IntegerVector child,
                   IntegerVector nTip) {
  if (parent.size() != child.size()) {
    throw std::length_error("Parent and child must be the same length");
  }
  const unsigned int n_tip = nTip[0],
                     n_internal = n_tip - 1,
                     n_edge = parent.size(),
                     all_node = n_internal + n_tip,
                     r_to_c = 1;
  if (n_edge != n_tip + n_tip - 2) {
    throw std::length_error("nEdge must == nTip + nTip - 2");
  }
  unsigned int smallest_below[MAX_NODE],
                    parent_of[MAX_NODE],
                    prime_id[MAX_NODE],
                    index[MAX_TIP];
  for (unsigned int i = 0; i != all_node; i++) {
    smallest_below[i] = i;
    prime_id[i] = i;
  }

  for (unsigned int i = 0; i != n_edge - 2; i += 2) {
    const unsigned int this_node = parent[i] - r_to_c,
      left_child = child[i] - r_to_c,
      right_child = child[i + 1] - r_to_c;
    smallest_below[this_node] = minimum(smallest_below[right_child],
                                        smallest_below[left_child]);
    prime_id[this_node] = maximum(smallest_below[left_child],
                              smallest_below[right_child]);
    parent_of[left_child] = parent_of[right_child] = this_node;

    for (unsigned int at = smallest_below[this_node]; at != this_node;
    at = parent_of[at]) {
      const unsigned int prime_candidate = prime_id[at];
      if (prime_candidate < prime_id[this_node]) {
        index[prime_id[this_node]] = prime_id[at] + (at < n_tip ? 0 : n_tip);
      }
    }
  }
  uint64_t num = 0U;
  uint64_t multiplier = 1U;
  for (unsigned int i = 3; i < n_tip; i++) {
    unsigned int insertion_edge = index[i];
    if (insertion_edge < n_tip) {
      --insertion_edge;
    } else {
      insertion_edge += i - (n_tip + 3);
    }
    num += insertion_edge * multiplier;
    multiplier *= (i + i - 3);
  }

  if (num >= MAX_INT) {
    NumericVector ret(2);
    ret[0] = num / MAX_INT;
    ret[1] = num % MAX_INT;
    return (ret);
  } else {
    NumericVector ret(1);
    ret[0] = num;
    return (ret);
  }
}
