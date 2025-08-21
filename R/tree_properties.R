#' Numeric index of each node in a tree
#' `NodeNumbers()` returns a sequence corresponding to the nodes in a tree
#' 
#' @template treeParam
#' @param tips Logical specifying whether to also include the indices of leaves.
#' @return `NodeNumbers()` returns an integer vector corresponding to the
#' indices of nodes within a tree.
#' @template MRS
#' @family tree properties
#' @family tree navigation
#' @export
NodeNumbers <- function(tree, tips = FALSE) {
  if (tips) {
    seq_len(NTip(tree) + tree[["Nnode"]])
  } else {
    NTip(tree) + seq_len(tree[["Nnode"]])
  }
}


#' Distance of each node from tree exterior
#'
#' `NodeDepth()` evaluates how "deep" each node is within a tree.
#'
#' For a rooted tree, the depth of a node is the minimum (if `shortest = TRUE`)
#' or maximum  (`shortest = FALSE`) number of edges that must be traversed,
#' moving away from the root, to reach a leaf.
#'
#' Unrooted trees are treated as if a root node occurs in the "middle" of the
#' tree, meaning the position that will minimise the maximum node depth.
#'
#'
#' @template xPhylo
#' @param shortest Logical specifying whether to calculate the length of the
#' shortest away-from-root path to a leaf.  If `FALSE`, the length of the
#' longest such route will be returned.
#' @param includeTips Logical specifying whether to include leaves
#' (each of depth zero) in return value.
#'
#' @return `NodeDepth()` returns an integer vector specifying the depth of
#' each external and internal node in `x`.
#'
#' @template MRS
#' @family tree navigation
#' @seealso [`ape::node.depth`] returns the number of tips descended from a
#' node.
#'
#' @examples
#' tree <- CollapseNode(BalancedTree(10), c(12:13, 19))
#' plot(tree)
#' nodelabels(NodeDepth(tree, includeTips = FALSE))
#'
#'
#' @export
NodeDepth <- function(x, shortest = FALSE, includeTips = TRUE) {
  UseMethod("NodeDepth")
}

#' @export
NodeDepth.list <- function(x, shortest = FALSE, includeTips = TRUE) {
  lapply(x, NodeDepth, shortest = shortest, includeTips = includeTips)
}

#' @export
NodeDepth.multiPhylo <- NodeDepth.list

#' @importFrom ape is.rooted
#' @export
NodeDepth.phylo <- function(x, shortest = FALSE, includeTips = TRUE) {
  if (is.rooted(x)) {
    .NodeDepth.rooted(x[["edge"]], shortest, includeTips)
  } else {
    NodeDepth(x[["edge"]], shortest, includeTips)
  }
}

#' @export
NodeDepth.matrix <- function(x, shortest = FALSE, includeTips = TRUE) {


  .NodeDepth.short <- function() {

    depths <- c(leaf0s, vapply(minVertex:nVertex, function(node)
      if (any(!is.na(leaf0s[child[parent == node]]))) 1L else NA_integer_
      , 0L))
    maxDepth <- 1L

    while(any(is.na(depths))) {
      for (node in rev(which(is.na(depths)))) {
        incident <- c(depths[child[parent == node]],
                      depths[parent[child == node]])
        na <- is.na(incident)
        nNa <- sum(na)
        if (nNa == 0L) {
          depths[node] <- min(incident) + 1L
        } else if (nNa == 1L) {
          aIncident <- incident[!na]
          if (all(aIncident <= maxDepth)) {
            depths[node] <- min(aIncident) + 1L
          }
        }
      }
      maxDepth <- maxDepth + 1L
    }

    #Return:
    depths
  }

  .NodeDepth.long <- function() {

    depths <- c(leaf0s, vapply(minVertex:nVertex, function(node)
      if (any(is.na(leaf0s[child[parent == node]]))) NA_integer_ else 1L
      , 0L))
    maxDepth <- 1L

    while(any(is.na(depths))) {
      for (node in rev(which(is.na(depths)))) {
        incident <- c(depths[child[parent == node]],
                      depths[parent[child == node]])
        na <- is.na(incident)
        nNa <- sum(na)
        if (nNa == 0L) {
          depths[node] <- sort(incident, decreasing = TRUE)[2] + 1L
        } else if (nNa == 1L) {
          aIncident <- incident[!na]
          if (all(aIncident <= maxDepth)) {
            depths[node] <- max(aIncident) + 1L
          }
        }
      }
      maxDepth <- maxDepth + 1L
    }

    #Return:
    depths
  }


  parent <- x[, 1]
  child <- x[, 2]
  minVertex <- min(parent)
  nVertex <- max(parent)

  nLeaf <- minVertex - 1L
  nNode <- nVertex - nLeaf
  leaf0s <- integer(nLeaf)

  depths <- if (shortest) .NodeDepth.short() else .NodeDepth.long()

  # Return:
  if (includeTips) depths else depths[minVertex:nVertex]

}

.NodeDepth.rooted <- function(x, shortest = FALSE, includeTips = TRUE) {

  parent <- x[, 1]
  child <- x[, 2]
  minVertex <- min(parent)
  nVertex <- max(parent)

  nLeaf <- minVertex - 1L
  nNode <- nVertex - nLeaf
  leaf0s <- integer(nLeaf)
  depths <- c(leaf0s, rep.int(NA_integer_, nNode))
  uncalculated <- is.na(depths)
  Func <- if (shortest) min else max

  while(any(uncalculated)) {
    depths[uncalculated] <- vapply(which(uncalculated), function(node) {
      Func(depths[child[parent == node]])
    }, 0L) + 1L
    uncalculated <- is.na(depths)
  }


  # Return:
  if (includeTips) depths else depths[minVertex:nVertex]

}

#' Number of edges incident to each node in a tree
#'
#' `NodeOrder()` calculates the order of each node: the number of edges
#' incident to it in a tree.
#' This value includes the root edge in rooted trees.
#'
#' @template xPhylo
#' @param includeAncestor Logical specifying whether to count edge leading to
#' ancestral node in calculation of order.
#' @param internalOnly Logical specifying whether to restrict to results
#' to internal nodes, i.e. to omit leaves. Irrelevant if
#' `includeAncestor = FALSE`.
#'
#' @return `NodeOrder()` returns an integer listing the order of each node;
#' entries are named with the number of each node.
#'
#' @examples
#' tree <- CollapseNode(BalancedTree(8), 12:15)
#' NodeOrder(tree)
#' plot(tree)
#' nodelabels(NodeOrder(tree, internalOnly = TRUE))
#'
#' @template MRS
#' @family tree navigation
#' @export
NodeOrder <- function(x, includeAncestor = TRUE, internalOnly = FALSE) {
  UseMethod("NodeOrder")
}


#' @export
NodeOrder.list <- function(x, includeAncestor = TRUE, internalOnly = FALSE) {
  lapply(x, NodeOrder, includeAncestor, internalOnly)
}

#' @export
NodeOrder.multiPhylo <- NodeOrder.list

#' @export
NodeOrder.phylo <- function(x, includeAncestor = TRUE, internalOnly = FALSE) {
  NodeOrder(x[["edge"]], includeAncestor, internalOnly)
}

#' @export
NodeOrder.matrix <- function(x, includeAncestor = TRUE, internalOnly = FALSE) {
  if (includeAncestor) {
    if (internalOnly) {
      tabulate(x)[-seq_len(min(x[, 1]) - 1L)]
    } else {
      tabulate(x)
    }
  } else {
    tabulate(x[, 1])[-seq_len(min(x[, 1]) - 1L)]
  }
}

#' Distance between edges
#'
#' Number of nodes that must be traversed to navigate from each edge to
#' each other edge within a tree
#'
#' @template treeParam
#'
#' @return `EdgeDistances()` returns a symmetrical matrix listing the number
#' of edges that must be traversed to travel from each numbered edge to each
#' other.
#' The two edges straddling the root of a rooted tree
#' are treated as a single edge.  Add a "root" tip using [`AddTip()`] if the
#' position of the root is significant.
#'
#' @examples
#'
#' tree <- BalancedTree(5)
#' plot(tree)
#' ape::edgelabels()
#'
#' EdgeDistances(tree)
#'
#' @family tree navigation
#' @importFrom fastmatch %fin%
#' @template MRS
#' @export
EdgeDistances <- function(tree) {
  edge <- tree[["edge"]]
  nEdge <- dim(edge)[1]
  edge <- RenumberTree(edge[, 1], edge[, 2])
  parent <- edge[, 1]
  child <- edge[, 2]
  ancs <- AllAncestors(parent, child)
  ret <- matrix(0L, ncol = nEdge, nrow = nEdge)
  for (i in seq_len(nEdge - 1L)) {
    ancI <- ancs[[child[i]]]
    nAncI <- length(ancI)
    for (j in seq.int(from = i + 1L, to = nEdge)) {
      ancJ <- ancs[[child[j]]]
      intersection <- intersect(ancI, ancJ)
      if (length(intersection) > 1L) {
        # On same side of root
        mrca <- max(intersection)
        if (child[i] %fin% ancJ || child[j] %fin% ancI) {
          addOne <- 0L
        } else {
          addOne <- 1L
        }
      } else {
        # On opposite sides of root
        mrca <- intersection
        addOne <- 0L
      }
      u <- union(ancI, ancJ)
      ret[i, j] <- sum(u > mrca, addOne)
    }
  }
  ret[lower.tri(ret)] <- t(ret)[lower.tri(ret)]

  origOrder <- match(tree[["edge"]][, 2], edge[, 2])

  # Return:
  ret[origOrder, origOrder]
}

#' Number of leaves in a phylogenetic tree
#'
#' `NTip()` extends [`ape::Ntip()`][ape::summary.phylo] to handle
#' objects of class `Splits` and `list`, and edge matrices
#' (equivalent to `tree$edge`).
#'
#' @param phy Object representing one or more phylogenetic trees.
#'
#' @return `NTip()` returns an integer specifying the number of tips in each
#' object in `phy`.
#'
#' @family tree properties
#' @export
NTip <- function(phy) UseMethod("NTip")

#' @rdname NTip
#' @export
NTip.default <- function(phy) attr(phy, "nTip")

#' @rdname NTip
#' @family Splits operations
#' @export
NTip.Splits <- NTip.default

#' @rdname NTip
#' @export
NTip.list <- function(phy) {
  vapply(phy, NTip, integer(1))
}

#' @rdname NTip
#' @export
NTip.phylo <- function(phy) {
  length(phy[["tip.label"]])
}

#' @rdname NTip
#' @export
NTip.multiPhylo <- function(phy) {
  vapply(phy, NTip.phylo, integer(1))
}

#' @rdname NTip
#' @export
NTip.phyDat <- function(phy) {
  length(phy)
}

#' @rdname NTip
#' @export
NTip.matrix <- function(phy) {
  if (is.numeric(phy)) {
    parent <- phy[, 1]
    child <- phy[, 2]

    # Return:
    max(child[!child %fin% parent])
  } else {
    NextMethod()
  }
}

#' Number of distinct splits
#'
#' `NSplits()` counts the unique bipartition splits in a tree or object.
#'
#' @param x A phylogenetic tree of class `phylo`; a list of such trees
#' (of class `list` or `multiPhylo`); a `Splits` object;
#' a vector of integers; or a character vector listing tips of a tree,
#' or a character of length one specifying a tree in Newick format.
#'
#' @return `NSplits()` returns an integer specifying the number of bipartitions in
#' the specified objects, or in a binary tree with `x` tips.
#'
#' @examples
#' NSplits(8L)
#' NSplits(PectinateTree(8))
#' NSplits(as.Splits(BalancedTree(8)))
#' @template MRS
#'
#' @family tree properties
#' @family Splits operations
#' @export
NSplits <- function(x) UseMethod("NSplits")

#' @rdname NSplits
#' @export
NPartitions <- NSplits

#' @rdname NSplits
#' @export
NSplits.phylo <- function(x) {
  # Use optimized C++ function that replaces collapse.singles bottleneck
  cpp_count_splits(x[["edge"]], length(x[["tip.label"]]))
}

#' @rdname NSplits
#' @export
NSplits.list <- function(x) vapply(x, NSplits, numeric(1L))

#' @rdname NSplits
#' @export
NSplits.multiPhylo <- NSplits.list

#' @rdname NSplits
#' @export
NSplits.Splits <- function(x) nrow(x)

#' @rdname NSplits
#' @export
NSplits.numeric <- function(x) pmax(0L, x - 3L)

#' @rdname NSplits
#' @export
NSplits.NULL <- function(x) NULL

#' @rdname NSplits
#' @export
NSplits.ClusterTable <- function(x) {
  # Root + Ingroup + All-leaves
  nrow(as.matrix(x)) - 3L
}

#' @rdname NSplits
#' @importFrom ape read.tree
#' @export
NSplits.character <- function(x) {
  lengthX <- length(x)
  if (lengthX == 1) {
    tree <- read.tree(text = x)
    if (!is.null(tree)) return(NSplits(tree))
  }
  NSplits(lengthX)
}

#' Maximum splits in an _n_-leaf tree
#'
#' `SplitsInBinaryTree()` is a convenience function to calculate the number of
#' splits in a fully-resolved (binary) tree with _n_ leaves.
#'
#' @param tree An object of a supported format that represents a tree or
#' set of trees, from which the number of leaves will be calculated.
#'
#' @return `SplitsInBinaryTree()` returns an integer vector detailing the number
#' of unique non-trivial splits in a binary tree with _n_ leaves.
#'
#' @examples
#' tree <- BalancedTree(8)
#' SplitsInBinaryTree(tree)
#' SplitsInBinaryTree(as.Splits(tree))
#' SplitsInBinaryTree(8)
#' SplitsInBinaryTree(list(tree, tree))
#' @template MRS
#'
#' @family tree properties
#' @family Splits operations
#' @export
SplitsInBinaryTree <- function(tree) UseMethod("SplitsInBinaryTree")

#' @rdname SplitsInBinaryTree
#' @export
SplitsInBinaryTree.list <- function(tree) vapply(tree, SplitsInBinaryTree, 0L)

#' @rdname SplitsInBinaryTree
#' @export
SplitsInBinaryTree.multiPhylo <- SplitsInBinaryTree.list

#' @rdname SplitsInBinaryTree
#' @export
SplitsInBinaryTree.numeric <- function(tree) as.integer(tree) - 3L

#' @rdname SplitsInBinaryTree
#' @export
SplitsInBinaryTree.NULL <- function(tree) NULL

#' @rdname SplitsInBinaryTree
#' @export
SplitsInBinaryTree.default <- function(tree) NTip(tree) - 3L

#' @rdname SplitsInBinaryTree
#' @export
SplitsInBinaryTree.Splits <- SplitsInBinaryTree.default

#' @rdname SplitsInBinaryTree
#' @export
SplitsInBinaryTree.phylo <- SplitsInBinaryTree.default

#' Is tree rooted?
#'
#' `TreeIsRooted()` is a fast alternative to `ape::is.rooted()`.
#'
#' @param tree A phylogenetic tree of class `phylo`.
#' @return `TreeIsRooted()` returns a logical specifying whether a root node is
#' resolved.
#'
#' @examples
#' TreeIsRooted(BalancedTree(6))
#' TreeIsRooted(UnrootTree(BalancedTree(6)))
#'
#' @family tree properties
#'
#' @template MRS
#' @export
TreeIsRooted <- function(tree) {
  edge <- tree[["edge"]]
  parent <- edge[, 1]
  sum(parent == min(parent)) < 3L
}

#' Which node is a tree's root?
#'
#' `RootNode()` identifies the root node of a (rooted or unrooted) phylogenetic
#' tree.
#' Unrooted trees are represented internally by a rooted tree with a polytomy
#' at the root.
#'
#' @param x A tree of class `phylo`, or its edge matrix; or a list or
#' `multiPhylo` object containing multiple trees.
#' @return `RootNode()` returns an integer denoting the root node for each tree.
#' Badly conformed trees trigger an error.
#' @template MRS
#'
#' @examples
#' RootNode(BalancedTree(8))
#' RootNode(UnrootTree(BalancedTree(8)))
#'
#'
#' @family tree navigation
#' @seealso
#'
#' Test whether a tree is rooted: [`TreeIsRooted()`]
#'
#' `phangorn::getRoot()`
#'
#' @export
RootNode <- function(x) UseMethod("RootNode")

#' @export
RootNode.phylo <- function(x) {
  edge <- x[["edge"]]
  edgeOrder <- attr(x, "order")
  if (!is.null(edgeOrder)) {
    if (edgeOrder == "postorder") {
      return(as.integer(edge[nrow(edge), 1L]))
    } else if (edgeOrder == "preorder") {
      return(as.integer(edge[1L]))
    }
  }
  RootNode(edge)
}

#' @export
RootNode.list <- function(x) {
  vapply(x, RootNode, 0L)
}

#' @export
RootNode.multiPhylo <- RootNode.list

#' @export
RootNode.numeric <- function(x) {
  parent <- x[, 1]
  child <- x[, 2]
  ret <- unique(parent[!parent %fin% child])
  if (length(ret) != 1) {
    warning("Root not unique: found ", paste(ret, collapse = ", "))
  }

  # Return:
  as.integer(ret)
}
