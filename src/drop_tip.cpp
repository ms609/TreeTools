#include <Rcpp.h>
using namespace Rcpp;

#include "../inst/include/TreeTools/types.h"
#include "../inst/include/TreeTools/renumber_tree.h"

#define NODE_COUNT(x) node_count[x - start_tip - 1]

#define RAISE_EDGE(i) if (n_deleted) ret((i) - n_deleted, _) = ret((i), _)

// [[Rcpp::export]]
IntegerMatrix drop_tip (const IntegerMatrix edge, const IntegerVector drop) {
  IntegerMatrix preorder = TreeTools::preorder_edges_and_nodes(edge(_, 0), edge(_, 1));

  const int32 root_node = preorder(0, 0),
    start_tip = root_node - 1,
    start_edge = preorder.nrow(),
    start_node = start_edge - start_tip
  ;

  int32
    still_to_drop = drop.length(),
    n_deleted = 0,
    ret_edges = start_edge - still_to_drop
  ;

  IntegerVector droppers = drop;
  IntegerMatrix ret(ret_edges, 2);

  std::unique_ptr<int32[]>
    node_count = std::make_unique<int32[]>(start_node)
  ;

  Rcout << "Starting on tree with " << start_tip << " tips, root = "
        << root_node <<", " << start_edge << " edges.\n";

  for (int32 i = start_edge; i--; ) {
    const int32 child = preorder(i, 1);
    bool delete_this = 0;
    if (child <= start_tip) {
      for (int32 j = still_to_drop; j--; ) {
        if (child == drop[j]) {
          droppers[j] = droppers[still_to_drop--];
          delete_this = 1;
          break;
        }
      }
    }
    if (delete_this) {
      n_deleted++;
      delete_this = 0;
    } else {
      ret(i - still_to_drop, _) = preorder(i, _);
      NODE_COUNT(preorder(i, 0);)++;
    }
  }

  do {
    n_deleted = 0;
    Rcout << ret_edges << " edges in ret\n";
    for (int32 i = 0; i != ret_edges; i++) { // postorder traversal
      Rcout << "   [" << i <<",]" << ret(i, 0) << "   " << ret(i, 1) << "\n";
      if (NODE_COUNT(ret(i, 1)) == 0) {
        Rcout << "Deleting dead end edge " << i <<"\n";
        NODE_COUNT(ret(i, 0))--;
        ++n_deleted;
      } else if (NODE_COUNT(ret(i, 1)) == 1) {
        Rcout << "Deleting singleton edge " << i <<"\n";
        NODE_COUNT(ret(i, 1))--;
        if (n_deleted) {
          ret(i - n_deleted, 0) = ret(i, 0);
          ret(i - n_deleted, 1) = ret(i + 1, 1);
        } else {
          ret(i, 1) = ret(i + 1, 1);
        }
        ++n_deleted;
      } else {
        RAISE_EDGE(i);
      }
    }
    ret_edges -= n_deleted;
  } while (n_deleted);

  return ret(Range(0, ret_edges), _);

/*

parent <- edge[, 1]
child <- edge[, 2]
external <- child <= nTip

# Drop tips:
  keep <- !child %in% which(drop)

# Drop dangling nodes:
    repeat {
  nonDanglers <- (child[keep] %in% parent[keep]) | external[keep]
  if (all(nonDanglers)) break
    keep[keep] <- nonDanglers
}

    parent <- parent[keep]
    child <- child[keep]

# Collapse singles:
    singletons <- tabulate(parent) == 1
    if (any(singletons)) {
      edgeBelowSingles <- rev(which(parent %in% which(singletons)))
      sortedSingles <- parent[edgeBelowSingles]
      edgeAboveSingles <- match(sortedSingles, child)

      for (i in seq_along(sortedSingles)) {
        child[edgeAboveSingles[i]] <- child[edgeBelowSingles[i]]
      }
      edge <- c(parent[-edgeBelowSingles], child[-edgeBelowSingles])
    } else {
      edge <- c(parent, child)
    }

    newNumbers <- integer(max(edge))
      uniqueInts <- unique(edge)
      newNumbers[uniqueInts] <- rank(uniqueInts)
      edge <- matrix(newNumbers[edge], ncol = 2L)
*/
}
