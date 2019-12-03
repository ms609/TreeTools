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

unsigned int n_options(const unsigned int a, const unsigned int b) {
  return n_shapes(a) * ((a == b) ? n_shapes(a + 1) / 2 : n_shapes(b));
}

unsigned int triangular_number(const unsigned int n) {
  return n * (n + 1) / 2;
}

unsigned int triangle_row (const unsigned int x) {
  // Solve x = (n)(n+1) / 2 for n
  return (sqrt((8 * x) + 1) - 1) / 2;
}

// [[Rcpp::export]]
IntegerVector n_rooted_shapes(IntegerVector nTip) {
  return IntegerVector::create(n_shapes(nTip[0]));
}

// Parent and child must be in postorder, with tree rooted.
// [[Rcpp::export]]
NumericVector edge_to_rooted_shape(IntegerVector parent, IntegerVector child,
                            IntegerVector nTip) {
  if (parent.size() != child.size()) {
    throw std::length_error("Parent and child must be the same length");
  }
  const unsigned int n_tip = nTip[0],
                     n_edge = parent.size(),
                     r_to_c = 1;
  if (n_edge != n_tip + n_tip - 2) {
    throw std::length_error("nEdge must == nTip + nTip - 2");
  }
  uint64_t tree_at[MAX_SHAPE_NODE] = {};
  unsigned int tips_below[MAX_SHAPE_NODE] = {};
  for (unsigned int i = 0; i != n_tip; i++) {
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

    const unsigned int small_tips = tips_below[small_child],
                       large_tips = tips_below[large_child];
    tips_below[this_node] = small_tips + large_tips;
    const unsigned int tips_here = tips_below[this_node];
    const unsigned int option_chosen = small_tips - 1;

    for (unsigned int unchosen_tips = 1; unchosen_tips <= option_chosen; unchosen_tips++) {
      tree_at[this_node] += n_options(unchosen_tips, tips_here - unchosen_tips);
    }

    if (small_tips == large_tips) {
      const unsigned int max_shape = n_shapes(large_tips),
        small_shape = tree_at[small_child],
        large_shape = tree_at[large_child];
      tree_at[this_node] += (triangular_number(max_shape) -
        triangular_number(max_shape - small_shape)) +
        large_shape - small_shape;
    } else {
      tree_at[this_node] += tree_at[large_child] +
      (tree_at[small_child] * n_shapes(large_tips));
    }
  }

  return NumericVector::create(tree_at[parent[n_edge - 1] - r_to_c]);
}

void fill_edges(unsigned int *parent, unsigned int *child,
                uint64_t n, const unsigned int n_tip,
                unsigned int *next_edge, unsigned int *next_tip,
                unsigned int *next_node) {

  const unsigned int this_node = (*next_node)++;

  for (unsigned int small_half = 1; ; small_half++) {
    const unsigned int large_half = n_tip - small_half;
    if (small_half == large_half) {
      parent[*next_edge] = this_node;
      if (small_half == 1) {
        child[(*next_edge)++] = (*next_tip)++;
        parent[*next_edge] = this_node;
        child[(*next_edge)++] = (*next_tip)++;
      } else {
        /*         LARGE TREE ID
         *          0  1   2  ...
         *  small 0 0  1   2  ... n
         *   tree 1   n+1 n+2 ... (2n - 1)
         *    id  2        2n ...
         */
        const unsigned int large_tree_options = n_shapes(large_half),
          total_options = triangular_number(large_tree_options),
          further_options = (total_options - n) - 1,
          small_tree = (large_tree_options - 1) - triangle_row(further_options),
          large_tree = small_tree + n - (total_options -
            triangular_number(triangle_row(further_options) + 1));
        child[(*next_edge)++] = *next_node;
        fill_edges(parent, child, small_tree, small_half,
                   next_edge, next_tip, next_node);

        parent[*next_edge] = this_node;
        child[(*next_edge)++] = *next_node;
        fill_edges(parent, child, large_tree, large_half,
                   next_edge, next_tip, next_node);
      }
      break;

    } else {
      const unsigned int small_trees = n_shapes(small_half),
        large_trees = n_shapes(large_half),
        options_here = small_trees * large_trees;

      if (n < options_here) {
        parent[*next_edge] = this_node;
        if (small_half == 1) {
          child[(*next_edge)++] = (*next_tip)++;
        } else {
          child[(*next_edge)++] = *next_node;
          fill_edges(parent, child, n / large_trees, small_half,
                     next_edge, next_tip, next_node);
        }

        parent[*next_edge] = this_node;
        if (large_half == 1) {
          child[(*next_edge)++] = (*next_tip)++;
        } else {
          child[(*next_edge)++] = *next_node;
          fill_edges(parent, child, n % large_trees, large_half,
                     next_edge, next_tip, next_node);
        }
        break;
      } else {
        n -= options_here;
      }
    }
  }

}

// [[Rcpp::export]]
IntegerMatrix rooted_shape_to_edge(NumericVector shape, IntegerVector nTip) {
  const unsigned int n_tip = nTip[0], n_edge = n_tip + n_tip - 2;
  uint64_t n = shape[0];
  unsigned int parent[MAX_SHAPE_NODE], child[MAX_SHAPE_NODE];
  unsigned int next_edge = 0, next_tip = 1, next_node = n_tip + 1;
  fill_edges(parent, child, n, n_tip, &next_edge, &next_tip, &next_node);

  IntegerMatrix ret (n_edge, 2);
  for (unsigned int i = 0; i < n_edge; i++) {
    ret(i, 0) = parent[i];
    ret(i, 1) = child[i];
  }
  return ret;
}
