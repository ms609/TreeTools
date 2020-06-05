#include <Rcpp.h>
#include "types.h"
using namespace Rcpp;
using namespace std;

const intx
  MAX_SHAPE_TIP = 55,
  MAX_SHAPE_NODE = MAX_SHAPE_TIP + MAX_SHAPE_TIP - 1
;

uint64_t n_shapes_cache[MAX_SHAPE_TIP + 1] = {};

uint64_t n_shapes(intx n_tips) {
  if (n_tips > 55) {
    throw std::length_error("64 bit integers cannot represent number of shapes"
                            "for > 55 tips");
  }
  if (!n_shapes_cache[n_tips]) {
    if (n_tips < 4) {
      n_shapes_cache[n_tips] = 1;
    } else {
      for (intx n_smaller = 1; n_smaller != ((n_tips + 1) / 2); n_smaller++) {
        const intx n_larger = n_tips - n_smaller;
        n_shapes_cache[n_tips] += (n_shapes(n_larger) * n_shapes(n_smaller));
      }
      if (n_tips % 2 == 0) {
        n_shapes_cache[n_tips] += ((n_shapes(n_tips / 2) *
          (n_shapes(n_tips / 2) + 1)) / 2);
      }
    }
  }
  return n_shapes_cache[n_tips];
}

uint64_t n_options(const intx a, const intx b) {
  return n_shapes(a) * ((a == b) ? n_shapes(a + 1) / 2 : n_shapes(b));
}

uint64_t triangular_number(const uint64_t n) {
  return n * (n + 1) / 2;
}

uint64_t triangle_row (const uint64_t x) {
  // Solve x = (n)(n+1) / 2 for n
  return uint64_t((sqrt((long double)(8 * x) + 1) - 1) / 2);
}

// [[Rcpp::export]]
IntegerVector n_rooted_shapes(IntegerVector nTip) {
  const uint64_t n = n_shapes(nTip[0]);
  if (n >= INT_MAX) {
    return IntegerVector{int(n / INT_MAX), int(n % INT_MAX)};
  } else {
    return IntegerVector{int(n)};
  }
}

// Parent and child must be in postorder, with tree rooted.
// [[Rcpp::export]]
IntegerVector edge_to_rooted_shape(IntegerVector parent, IntegerVector child,
                            IntegerVector nTip) {
  if (parent.length() != child.length()) {
    throw std::length_error("Parent and child must be the same length");
  }
  const intx
    n_tip = nTip[0],
    n_edge = parent.length(),
    r_to_c = 1
  ;
  if (n_tip > 55) {
    // Limit of 64-bit integers
    throw std::range_error("Cannot calculate shape with > 55 leaves");
  }
  if (n_edge != n_tip + n_tip - 2) {
    throw std::length_error("nEdge must == nTip + nTip - 2: is tree binary?");
  }
  uint64_t tree_at[MAX_SHAPE_NODE] = {};
  intx tips_below[MAX_SHAPE_NODE] = {};
  for (intx i = 0; i != n_tip; i++) {
    tree_at[i] = 0;
    tips_below[i] = 1;
  }

  for (intx i = 0; i != n_edge; i += 2) {
    const intx
      this_node = parent[i] - r_to_c,
      left_child = child[i] - r_to_c,
      right_child = child[i + 1] - r_to_c
    ;
    intx large_child, small_child;

    if ((tips_below[left_child] > tips_below[right_child]) ||
      ((tips_below[left_child] == tips_below[right_child]) &&
      tree_at[left_child] > tree_at[right_child])) {
      small_child = right_child;
      large_child = left_child;
    } else {
      small_child = left_child;
      large_child = right_child;
    }

    const intx
      small_tips = tips_below[small_child],
      large_tips = tips_below[large_child]
    ;
    tips_below[this_node] = small_tips + large_tips;
    const intx
      tips_here = tips_below[this_node],
      option_chosen = small_tips - 1
    ;

    for (intx unchosen_tips = 1;
         unchosen_tips <= option_chosen;
         unchosen_tips++) {
      tree_at[this_node] += n_options(unchosen_tips, tips_here - unchosen_tips);
    }

    if (small_tips == large_tips) {
      const uint64_t
        max_shape = n_shapes(large_tips),
        small_shape = tree_at[small_child],
        large_shape = tree_at[large_child]
      ;
      tree_at[this_node] += (triangular_number(max_shape) -
        triangular_number(max_shape - small_shape)) +
        large_shape - small_shape;
    } else {
      tree_at[this_node] += tree_at[large_child] +
      (tree_at[small_child] * n_shapes(large_tips));
    }
  }

  const uint64_t ret = tree_at[parent[n_edge - 1] - r_to_c];
  if (ret >= INT_MAX) {
    return IntegerVector{int(ret / INT_MAX), int(ret % INT_MAX)};
  } else {
    return IntegerVector{int(ret)};
  }
}

void fill_edges(intx *parent, intx *child, uint64_t n, const intx n_tip,
                intx *next_edge, intx *next_tip, intx *next_node) {

  const intx this_node = (*next_node)++;

  for (intx small_half = 1; ; small_half++) {

    const intx large_half = n_tip - small_half;

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
        const uint64_t
          large_tree_options = n_shapes(large_half),
          total_options = triangular_number(large_tree_options),
          further_options = (total_options - n) - 1,
          small_tree = (large_tree_options - 1) - triangle_row(further_options),
          large_tree = small_tree + n - (total_options -
            triangular_number(triangle_row(further_options) + 1))
        ;
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
      const uint64_t
        small_trees = n_shapes(small_half),
        large_trees = n_shapes(large_half),
        options_here = small_trees * large_trees
      ;

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
  const intx
    n_tip = nTip[0],
    n_edge = n_tip + n_tip - 2;
  intx
    parent[MAX_SHAPE_NODE],
    child[MAX_SHAPE_NODE],
    next_edge = 0,
    next_tip = 1,
    next_node = n_tip + 1;
  uint64_t n = shape[0];
  fill_edges(parent, child, n, n_tip, &next_edge, &next_tip, &next_node);

  IntegerMatrix ret (n_edge, 2);
  for (intx i = 0; i != n_edge; i++) {
    ret(i, 0) = parent[i];
    ret(i, 1) = child[i];
  }
  return ret;
}
