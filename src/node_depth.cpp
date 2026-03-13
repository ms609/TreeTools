#include <Rcpp/Lightest>
#include <algorithm> /* for min_element, max_element, sort */
#include <vector>
using namespace Rcpp;

// Compute node depth for unrooted trees in O(n) via two-pass traversal.
// depth[v] = second-largest neighbour distance (shortest=false)
//         or smallest neighbour distance (shortest=true)
// Tips have depth 0.
// Returns depths for all vertices 1..max_node.
// [[Rcpp::export]]
IntegerVector node_depth_unrooted(
    const IntegerVector parent,
    const IntegerVector child,
    const IntegerVector postorder, // 1-based edge indices in postorder
    const bool shortest
) {
  const int n_edge = parent.size();
  const int root = *std::min_element(parent.begin(), parent.end());
  const int n_tip = root - 1;
  const int max_node = *std::max_element(parent.begin(), parent.end());

  // subtree_h[v]: max (or min) depth below v in the rooted representation
  // Tips: 0.  Internal: func(children's subtree_h + 1)
  std::vector<int> subtree_h(max_node + 1, 0);

  // For longest path: store top-2 subtree heights per node
  // For shortest path: store bottom-2
  // best1[v] >= best2[v] (longest) or best1[v] <= best2[v] (shortest)
  // best1_child[v] = which child contributed best1
  const int INF = 1000000;
  const int NEG_INF = -1;
  std::vector<int> best1(max_node + 1, shortest ? INF : NEG_INF);
  std::vector<int> best2(max_node + 1, shortest ? INF : NEG_INF);
  std::vector<int> best1_child(max_node + 1, -1);

  // Postorder: compute subtree_h and best1/best2
  for (int i = 0; i < n_edge; ++i) {
    const int e = postorder[i] - 1; // 0-based
    const int p = parent[e];
    const int c = child[e];
    const int ch = subtree_h[c] + 1; // height through this child

    if (shortest) {
      if (ch <= best1[p]) {
        best2[p] = best1[p];
        best1[p] = ch;
        best1_child[p] = c;
      } else if (ch < best2[p]) {
        best2[p] = ch;
      }
      subtree_h[p] = best1[p];
    } else {
      if (ch >= best1[p]) {
        best2[p] = best1[p];
        best1[p] = ch;
        best1_child[p] = c;
      } else if (ch > best2[p]) {
        best2[p] = ch;
      }
      subtree_h[p] = best1[p];
    }
  }

  // from_above[v]: longest (or shortest) path from v going through parent
  std::vector<int> from_above(max_node + 1, NEG_INF);
  // Root has no parent; from_above[root] stays at NEG_INF

  // Preorder: compute from_above (reverse postorder)
  for (int i = n_edge - 1; i >= 0; --i) {
    const int e = postorder[i] - 1;
    const int p = parent[e];
    const int c = child[e];

    // Best alternative at p excluding path through c:
    // If c contributed best1[p], use best2[p]; otherwise use best1[p]
    int best_sibling;
    if (best1_child[p] == c) {
      best_sibling = best2[p];
    } else {
      best_sibling = best1[p];
    }

    if (shortest) {
      // from_above[c] = 1 + min(from_above[p], best_sibling)
      if (from_above[p] == NEG_INF) {
        // p is root (or from_above not set): only siblings matter
        from_above[c] = (best_sibling < INF) ? 1 + best_sibling : INF;
      } else {
        from_above[c] = 1 + std::min(from_above[p], best_sibling);
      }
    } else {
      if (from_above[p] == NEG_INF) {
        from_above[c] = (best_sibling > NEG_INF) ? 1 + best_sibling : NEG_INF;
      } else {
        from_above[c] = 1 + std::max(from_above[p], best_sibling);
      }
    }
  }

  // Combine: depth[v] for internal nodes
  IntegerVector depths(max_node, 0); // 1-indexed via [v-1]; tips stay 0

  for (int v = n_tip + 1; v <= max_node; ++v) {
    // Collect all neighbour distances
    // Children contribute subtree_h[child] + 1 (already in best1/best2)
    // Parent direction contributes from_above[v]

    if (shortest) {
      int d = best1[v]; // smallest child path
      if (from_above[v] != NEG_INF && from_above[v] < d) {
        d = from_above[v];
      }
      depths[v - 1] = d;
    } else {
      // Need second-largest among all neighbour distances
      // Neighbours: best1[v], best2[v], from_above[v]
      // (best1 >= best2 for longest)
      // Sort descending and take [1] (second largest)
      int a = best1[v];
      int b = best2[v];
      int c_val = from_above[v];
      // Ensure we handle cases where from_above is undefined (root)
      if (c_val == NEG_INF) {
        // Root node: only children matter; second largest = best2
        depths[v - 1] = (b > NEG_INF) ? b : a;
      } else {
        // Three candidates: a >= b, c_val unknown
        if (c_val >= a) {
          depths[v - 1] = a; // c_val >= a >= b, second = a
        } else if (c_val >= b) {
          depths[v - 1] = c_val; // a > c_val >= b, second = c_val
        } else {
          depths[v - 1] = b; // a >= b > c_val, second = b
        }
      }
    }
  }

  return depths;
}
