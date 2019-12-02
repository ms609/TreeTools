#include <stdint.h>
#include <Rcpp.h>
using namespace Rcpp;


const unsigned int MAX_SHAPE_TIP = 200,
  MAX_SHAPE_NODE = MAX_SHAPE_TIP + MAX_SHAPE_TIP - 1;

unsigned int n_shapes_cache[MAX_SHAPE_TIP + 1] = {};

unsigned int n_shapes (unsigned int n) {
  if (!n_shapes_cache[n]) {
    if (n < 4) {
      n_shapes_cache[n] = 1;
    } else {
      for (unsigned int n_smaller = 1; n_smaller < ((n + 1) / 2); n_smaller++) {
        const unsigned int n_larger = n - n_smaller;
        n_shapes_cache[n] += (n_shapes(n_larger) * n_shapes(n_smaller));
      }
      if (n % 2 == 0) {
        n_shapes_cache[n] += ((n_shapes(n/2) * (n_shapes(n/2) + 1)) / 2);
      }
    }
  }
  return n_shapes_cache[n];
}

// Parent and child must be in postorder, with tree rooted on tip 1.
// [[Rcpp::export]]
NumericVector edge_to_shape(IntegerVector parent, IntegerVector child,
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
  uint64_t tree_at[MAX_SHAPE_NODE] = {};
  unsigned int tips_below[MAX_SHAPE_NODE] = {};
  for (unsigned int i = 0; i != all_node; i++) {
    tree_at[i] = 0;
    tips_below[i] = 1;
  }

  for (unsigned int i = 0; i != n_edge; i += 2) {
    const unsigned int this_node = parent[i] - r_to_c,
      left_child = child[i] - r_to_c,
      right_child = child[i + 1] - r_to_c;
    unsigned int large_child, small_child;


    if ((tips_below[left_child] > tips_below[right_child]) ||
      ((tips_below[left_child] == tips_below[right_child]) &&
      tree_at[left_child] > tree_at[right_child])) {
      small_child = right_child;
      large_child = left_child;
    } else {
      small_child = left_child;
      large_child = right_child;
    }

    tips_below[this_node] = tips_below[small_child] + tips_below[large_child];
    const unsigned int options_here = tips_below[this_node] / 2;
    const unsigned int option_chosen = tips_below[small_child] - 1;
    tree_at[this_node] = option_chosen +
      (tree_at[small_child] * options_here) +
      (tree_at[large_child] * options_here * n_shapes(tips_below[small_child]));
  }
  return NumericVector::create(tree_at[parent[n_edge - 1] - r_to_c]);
}
