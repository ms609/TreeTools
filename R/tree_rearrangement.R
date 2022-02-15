#' Root or unroot a phylogenetic tree
#'
#' `RootTree()` roots a tree on the smallest clade containing the specified
#' tips;
#' `RootOnNode()` roots a tree on a specified internal node;
#' `UnrootTree()` collapses a root node, without the undefined behaviour
#' encountered when using \code{\link[ape:root]{ape::unroot}()} on trees in
#' preorder.
#'
#' @template tree(s)Param
#' @param outgroupTips Vector of type character, integer or logical, specifying
#' the names or indices of the tips to include in the outgroup.  If
#' `outgroupTips` is a of type character, and a tree contains multiple tips
#' with a matching label, the first will be used.
#'
#' @return `RootTree()` returns a tree of class `phylo`, rooted on the smallest
#' clade that contains the specified tips, with edges and nodes numbered in
#' preorder.
#'
#' @examples
#' tree <- PectinateTree(8)
#' plot(tree)
#' ape::nodelabels()
#'
#' plot(RootTree(tree, c('t6', 't7')))
#'
#' plot(RootOnNode(tree, 12))
#' plot(RootOnNode(tree, 2))
#'
#' @seealso
#' - [`ape::root()`]
#' - [`EnforceOutgroup()`]
#'
#' @family tree manipulation
#'
#' @template MRS
#' @export
RootTree <- function(tree, outgroupTips) {
  if (missing(outgroupTips)) return(tree)
  if (is.null(outgroupTips) ||
      length(outgroupTips) == 0) {
    return(tree)
  }
  UseMethod('RootTree')
}

#' @export
RootTree.phylo <- function(tree, outgroupTips) {
  if (missing(outgroupTips) ||
      is.null(outgroupTips) ||
      length(outgroupTips) == 0) {
    return(tree)
  }
  tipLabels <- tree[["tip.label"]]
  if (is.character(outgroupTips)) {
    chosenTips <- match(outgroupTips, tipLabels)
    if (any(is.na(chosenTips))) {
      stop("Outgroup tips [",
           paste(outgroupTips[is.na(chosenTips)], collapse = ', '),
           "] not found in tree's tip labels.")
    }
    outgroupTips <- chosenTips
  } else if (is.logical(outgroupTips)) {
    outgroupTips <- which(outgroupTips)
  }
  tree <- Preorder(tree)
  nTip <- NTip(tree)
  if (length(outgroupTips) == 0) {
    stop("No outgroup tips selected")
  } else if (length(outgroupTips) == 1L) {
    outgroup <- outgroupTips
  } else if (length(outgroupTips) == nTip - 1L) {
    outgroup <- setdiff(seq_len(nTip), outgroupTips)
  } else {
    ancestryTable <- .AncestorTable(tree, outgroupTips)
    lineage <- ancestryTable[1, ]
    lca <- max(lineage[ancestryTable[2, ] == length(outgroupTips)])
    nTip <- length(tipLabels)
    rootNode <- nTip + 1L
    if (lca == rootNode) {
      if (tree[["Nnode"]] > 2L) {
        lca <- lineage[lineage - c(lineage[-1], 0) != -1][1] + 1L
        if (lca > rootNode + nTip - 2L) {
          return(tree)
        }
      } else {
        return(tree)
      }
    }
    outgroup <- lca
  }
  if (outgroup > nTip + tree[["Nnode"]]) {
    return(tree)
  }

  root_on_node(tree, outgroup)
}

.AllAncestors <- function(edge) {
  edge <- Preorder(edge)
  parents <- edge[, 1]
  child <- edge[, 2]
  nEdge <- length(child)
  res <- vector("list", max(parents))
  for (i in seq_len(nEdge)) {
    pa <- parents[i]
    res[[child[i]]] <- c(pa, res[[pa]])
  }

  # Return:
  res
}

.AncestorTable <- function(tree, outgroupTips) {
  edge <- tree[["edge"]]
  parent <- edge[, 1]
  child <- edge[, 2]
  nVert <- max(parent)
  parentOf <- rep_len(NA_integer_, nVert)
  parentOf[child] <- parent
  
  counts <- integer(nVert)
  i <- outgroupTips
  while (length(i)) {
    i <- parentOf[i]
    i <- i[!is.na(i)]
    for (j in i) {
      counts[j] <- counts[j] + 1L
    }
  }
  
  counted <- counts > 0
  
  # Return:
  rbind(node = which(counted),
        count = counts[counted])
}

#' @export
RootTree.matrix <- function(tree, outgroupTips) {
  tree <- Preorder(tree)
  if (missing(outgroupTips) ||
      is.null(outgroupTips) ||
      length(outgroupTips) == 0L) {
    return(tree)
  }
  rootNode <- tree[1]
  nNode <- max(tree[, 1]) - rootNode + 1L
  if (length(outgroupTips) == 1L) {
    outgroup <- outgroupTips
  } else {
    ancestry <- unlist(.AllAncestors(tree)[outgroupTips])
    ancestryTable <- tabulate(ancestry)
    lca <- max(which(ancestryTable == length(outgroupTips)))
    
    if (lca == rootNode) {
      if (nNode > 2L) {
        lineage <- which(as.logical(ancestryTable))
        lca <- lineage[lineage - c(lineage[-1], 0) != -1][1] + 1L
      } else {
        return(tree)
      }
    }
    outgroup <- lca
  }

  root_on_node(list(edge = tree, Nnode = nNode), outgroup)[["edge"]]
}

#' @export
RootTree.list <- function(tree, outgroupTips) {
  lapply(tree, RootTree, outgroupTips)
}

#' @export
RootTree.multiPhylo <- function(tree, outgroupTips) {
  tree[] <- RootTree.list(tree, outgroupTips)
  tree
}

#' @export
RootTree.NULL <- function(tree, outgroupTips) NULL

#' @rdname RootTree
#' @param node integer specifying node (internal or tip) to set as the root.
#' @param resolveRoot logical specifying whether to resolve the root node.
#'
#' @return `RootOnNode()` returns a tree of class `phylo`, rooted on the
#' requested `node` and ordered in [`Preorder`].
#'
#' @export
RootOnNode <- function(tree, node, resolveRoot = FALSE) {
  UseMethod('RootOnNode', tree)
}

#' @importFrom fastmatch %fin%
#' @export
RootOnNode.phylo <- function(tree, node, resolveRoot = FALSE) {
  edge <- tree[["edge"]]
  parent <- edge[, 1]
  child <- edge[, 2]
  nodeParentEdge <- child == node
  rootEdges <- !parent %fin% child
  rootNode <- parent[rootEdges][1]
  rooted <- sum(rootEdges) == 2L

  if (any(nodeParentEdge)) {
    if (rooted) {
      # Do before editing parent:
      ancestorEdges <- EdgeAncestry(which(nodeParentEdge), parent, child,
                                    stopAt = rootEdges)

      rootChildren <- child[rootEdges]
      if (node %fin% rootChildren) {
        if (!resolveRoot) {
          if (node < NTip(tree)) {
            node <- rootChildren[rootChildren != node]
          }
          deletedEdge  <- child == node
          parent[parent == node] <- parent[deletedEdge]
          parent <- parent[!deletedEdge]
          child <- child[!deletedEdge]
          parent[parent > node] <- parent[parent > node] - 1L
          child[child > node] <- child[child > node] - 1L
          tree[["Nnode"]] <- tree[["Nnode"]] - 1L
        }
        inverters <- logical(length(parent))
      } else {
        spareRoot <- rootChildren == min(rootChildren) # Hit tip if present
        parent[rootEdges][!spareRoot] <- rootChildren[spareRoot]

        inverters <- ancestorEdges | rootEdges
        if (!ancestorEdges[rootEdges][!spareRoot]) {
          inverters[which(rootEdges)[!spareRoot]] <- FALSE
        }
        if (resolveRoot) {
          inverters <- inverters | nodeParentEdge
          parent[rootEdges][spareRoot] <- node
          child[rootEdges][spareRoot] <- rootNode
          child[nodeParentEdge] <- rootNode
        } else {
          if (node > rootNode) inverters <- inverters | nodeParentEdge
          deletedEdge <- which(rootEdges)[spareRoot]

          parent <- parent[-deletedEdge]
          child <- child[-deletedEdge]
          inverters <- inverters[-deletedEdge]

          parent[parent > rootNode] <- parent[parent > rootNode] - 1L
          child[child > rootNode] <- child[child > rootNode] - 1L

          tree[["Nnode"]] <- tree[["Nnode"]] - 1L
        }
      }
    } else {
      if (resolveRoot) {
        inverters <- c(EdgeAncestry(which(nodeParentEdge), parent, child), FALSE)
        newNode <- max(parent) + 1L
        parent <- c(parent, newNode)
        child <- c(child, parent[nodeParentEdge])
        parent[nodeParentEdge] <- newNode
        tree[["Nnode"]] <- 1L + tree[["Nnode"]]
      } else {
        inverters <- EdgeAncestry(which(nodeParentEdge), parent, child,
                                  stopAt = rootEdges)
      }
    }
    tree[["edge"]] <- RenumberTree(ifelse(inverters, child, parent),
                                   ifelse(inverters, parent, child))
    attr(tree, 'order') <- 'preorder'
    tree
  } else {
    # Root position is already correct
    if (rooted && !resolveRoot) {
      UnrootTree(tree)
    } else if (!rooted && resolveRoot) {
      RootOnNode(tree, max(child[rootEdges]), resolveRoot = TRUE)
    } else {
      Preorder(tree)
    }
  }
}

#' @export
RootOnNode.list <- function(tree, node, resolveRoot = FALSE) {
  lapply(tree, RootOnNode, node, resolveRoot)
}

#' @export
RootOnNode.multiPhylo <- function(tree, node, resolveRoot = FALSE) {
  tree[] <- RootOnNode.list(tree, node, resolveRoot)
  tree
}

#' @export
RootOnNode.NULL <- function(tree, node, resolveRoot = FALSE) NULL

#' @rdname RootTree
#' @return `UnrootTree()` returns `tree`, in preorder,
#' having collapsed the first child of the root node in each tree.
#' @export
UnrootTree <- function(tree) UseMethod('UnrootTree')

#' @export
UnrootTree.phylo <- function(tree) {
  tree <- Preorder(tree)
  edge <- tree[["edge"]]
  if (dim(edge)[1] < 3) {
    return(tree)
  }

  parent <- edge[, 1]
  rootNode <- parent[1]
  rootEdge2 <- parent[-1] == rootNode
  if (sum(rootEdge2) > 1L) {
    return(tree)
  }

  if (edge[1, 2] < rootNode) {
    deleted <- 1L + which(rootEdge2)
    weightTo <- 1L
  } else {
    deleted <- 1L
    weightTo <- 1L + which(rootEdge2)
  }
  renumber <- edge > rootNode
  edge[renumber] <- edge[renumber] - 1L
  tree[["edge"]] <- edge[-deleted, ]
  tree[["Nnode"]] <- tree[["Nnode"]] - 1L
  weight <- tree[["edge.length"]]
  if (!is.null(weight)) {
    weight[weightTo] <- weight[weightTo] + weight[deleted]
    tree[["edge.length"]] <- weight[-deleted]
  }

  # Return:
  tree
}

#' @export
UnrootTree.list <- function(tree) lapply(tree, UnrootTree)

#' @export
UnrootTree.multiPhylo <- function(tree) {
  tree[] <- UnrootTree.list(tree)
  tree
}

#' @export
UnrootTree.NULL <- function(tree) NULL

#' Collapse nodes on a phylogenetic tree
#'
#' Collapses specified nodes or edges on a phylogenetic tree, resulting in
#' polytomies.
#'
#' @template treeParam
#' @param nodes,edges Integer vector specifying the nodes or edges in the tree
#'  to be dropped.
#' (Use \code{\link[ape]{nodelabels}()} or
#' \code{\link[ape:nodelabels]{edgelabels}()}
#' to view numbers on a plotted tree.)
#'
#' @return `CollapseNode()` and `CollapseEdge()` return a tree of class `phylo`,
#' corresponding to `tree` with the specified nodes or edges collapsed.
#' The length of each dropped edge will (naively) be added to each descendant
#' edge.
#'
#' @examples
#' oldPar <- par(mfrow = c(3, 1), mar = rep(0.5, 4))
#'
#' tree <- as.phylo(898, 7)
#' tree$edge.length <- 11:22
#' plot(tree)
#' nodelabels()
#' edgelabels()
#' edgelabels(round(tree$edge.length, 2),
#'            cex = 0.6, frame = 'n', adj = c(1, -1))
#'
#' # Collapse by node number
#' newTree <- CollapseNode(tree, c(12, 13))
#' plot(newTree)
#' nodelabels()
#' edgelabels(round(newTree$edge.length, 2),
#'            cex = 0.6, frame = 'n', adj = c(1, -1))
#'
#' # Collapse by edge number
#' newTree <- CollapseEdge(tree, c(2, 4))
#' plot(newTree)
#'
#' par(oldPar)
#'
#' @family tree manipulation
#' @author  Martin R. Smith
#' @export
CollapseNode <- function(tree, nodes) UseMethod('CollapseNode')

#' @rdname CollapseNode
#' @importFrom fastmatch %fin%
#' @export
CollapseNode.phylo <- function(tree, nodes) {
  if (length(nodes) == 0) return(tree)

  edge <- tree[["edge"]]
  lengths <- tree[["edge.length"]]
  hasLengths <- !is.null(lengths)
  parent <- edge[, 1]
  child <- edge[, 2]
  root <- RootNode(edge)
  nTip <- NTip(tree)
  maxNode <- max(parent)
  edgeBelow <- order(child, method = 'radix') # a little faster than 'auto'
  tips <- seq_len(nTip)
  preRoot <- seq_len(root - 1L)
  edgeBelow <- c(edgeBelow[preRoot], NA, edgeBelow[-preRoot])
  nodes <- unique(nodes)

  if (!all(nodes %fin% (nTip + 1L):maxNode)) {
    stop("nodes must be integers between ", nTip + 1L, " and ", maxNode)
  }
  if (any(nodes == root)) {
    warning("Cannot collapse root node")
    nodes <- nodes[nodes != root]
    if (length(nodes) == 0L) return(tree)
  }

  keptEdges <- -edgeBelow[nodes]
  depths <- NodeDepth(edge)[nodes]

  for (node in nodes[order(depths)]) {
    newParent <- parent[edgeBelow[node]]
    if (hasLengths) lengths[parent == node] <- lengths[parent == node] + lengths[child == node]
    parent[parent == node] <- newParent
  }

  newNumber <- c(seq_len(nTip), nTip + cumsum((nTip + 1L):maxNode %fin% parent))

  tree[["edge"]] <- cbind(newNumber[parent[keptEdges]],
                          newNumber[child[keptEdges]])
  tree[["edge.length"]] <- lengths[keptEdges]
  tree[["Nnode"]] <- tree[["Nnode"]] - length(nodes)

  # TODO Renumber nodes sequentially
  # TODO Re-write this in C++.
  tree
}

#' @rdname CollapseNode
#' @export
CollapseEdge <- function(tree, edges) {
  nodesToCollapse <- tree[["edge"]][edges, 2]
  if (any(nodesToCollapse < NTip(tree))) {
    stop("Cannot collapse external edges: ",
         paste(edges[nodesToCollapse <= NTip(tree)], collapse = ', '))
  }
  CollapseNode(tree, nodesToCollapse)
}

#' Generate binary tree by collapsing polytomies
#'
#' `MakeTreeBinary()` resolves, at random, all polytomies in a tree or set of
#' trees, such that all trees compatible with the input topology are drawn
#' with equal probability.
#'
#' @seealso Since ape v5.5, this functionality is available through
#' [`ape::multi2di()`]; previous versions of 'ape' did not return topologies
#' in equal frequencies.
#'
#' @return `MakeTreeBinary()` returns a rooted binary tree of class `phylo`,
#' corresponding to tree uniformly selected from all those compatible with
#' the input tree topologies.
#'
#' @examples
#' MakeTreeBinary(CollapseNode(PectinateTree(7), c(9, 11, 13)))
#' UnrootTree(MakeTreeBinary(StarTree(5)))
#' @template MRS
#' @template treeParam
#' @family tree manipulation
#' @export
MakeTreeBinary <- function(tree) {
  UseMethod('MakeTreeBinary')
}

#' @export
MakeTreeBinary.phylo <- function(tree) {
  tree <- Preorder(tree)
  degree <- NodeOrder(tree, internalOnly = TRUE)
  degree[1] <- degree[1] + 1L # Root node
  polytomies <- degree > 3L
  if (!any(polytomies)) return(tree)
  edge <- tree[["edge"]]

  nTip <- edge[1] - 1L
  polytomyN <- which(polytomies) + nTip
  degree <- degree[polytomies]
  for (i in seq_len(sum(polytomies))) {
    n <- polytomyN[i]
    nKids <- degree[i] - 1L
    newParent <- .RandomParent(nKids + 1L) # Tip 1 is the 'root'
    newEdges <- RenumberEdges(newParent, seq_len(nKids + nKids))

    nNewNodes <- nKids - 2L

    childEdges <- edge[, 1] == n

    keep <- edge[!childEdges, ]
    increase <- keep > n
    keep[increase] <- keep[increase] + nNewNodes

    children <- edge[childEdges, 2L]
    increase <- children > n
    children[increase] <- children[increase] + nNewNodes

    polytomyN <- polytomyN + nNewNodes

    newEdges2 <- newEdges[[2]][c(-1, -2)] - 1L
    decrease <- newEdges2 > nKids
    newEdges2[decrease] <- newEdges2[decrease] - 2L

    add <- cbind(newEdges[[1]][-(1:2)] + n - nKids - 3L,
                 c(children, n + seq_len(nNewNodes))[newEdges2])

    edge <- rbind(keep, add)
  }
  tree[["edge"]] <- edge
  tree[["Nnode"]] <- nTip - 1L
  tree
}

#' @export
MakeTreeBinary.list <- function(tree) lapply(tree, MakeTreeBinary)

#' @export
MakeTreeBinary.multiPhylo <- function(tree) {
  structure(MakeTreeBinary.list(tree), class = 'multiPhylo')
}


#' Leaf label interchange
#'
#' `LeafLabelInterchange()` exchanges the position of leaves within a tree.
#'
#' Modifies a tree by switching the positions of _n_ leaves.  To avoid
#' later swaps undoing earlier exchanges, all _n_ leaves are guaranteed
#' to change position.  Note, however, that no attempt is made to avoid
#' swapping equivalent leaves, for example, a pair that are each others'
#' closest relatives.  As such, the relationships within a tree are not
#' guaranteed to be changed.
#'
#' @template treeParam
#' @param n Integer specifying number of leaves whose positions should be
#' exchanged.
#'
#' @return `LeafLabelInterchange()` returns a tree of class `phylo` on which
#' the position of `n` leaves have been exchanged.
#' The tree's internal topology will not change.
#'
#' @examples
#' tree <- PectinateTree(8)
#' plot(LeafLabelInterchange(tree, 3L))
#'
#' @template MRS
#' @family tree manipulation
#' @export
LeafLabelInterchange <- function(tree, n = 2L) {

  if (n < 2L) return(tree)
  tipLabel <- tree[["tip.label"]]
  nTip <- length(tipLabel)
  if (n > nTip) {
    stop("`n` cannot exceed number of tips in tree (", nTip, ")")
  }

  from <- sample.int(nTip, n)
  # We need to ensure that a tip's new position is different from its old one.
  # A simple way would be to shift each moved tip to the next one in turn:
  # 1, 2, 3, 4 -> 2, 3, 4, 1
  # But this cyclicity is non-random: we may wish to introduce some smaller
  # cycles, e.g.
  # 1, 2, 3, 4 -> 2, 1, 4, 3.

  cycles <- integer(n / 2L)
  remaining <- n
  i <- 1L
  while (remaining > 0L) {
    if (remaining == 2L) {
      cycles[i] <- 2L
      break
    }
    thisCycle <- SampleOne(1L + seq_len(remaining - 2L))
    if (thisCycle == remaining - 1L) thisCycle <- remaining
    remaining <- remaining - thisCycle
    cycles[i] <- thisCycle
    i <- i + 1L
  }
  cycles <- cycles[cycles > 0L]
  nCycles <- length(cycles)
  start <- cumsum(c(0L, cycles[-nCycles]))
  to <- unlist(lapply(seq_along(cycles), function(i) {
     c(seq_len(cycles[i] - 1L) + 1L, 1L) + start[i]
  }))

  tree[["tip.label"]][from] <- tipLabel[from[to]]
  tree
}

