#' Add a tip to a phylogenetic tree
#'
#' `AddTip()` adds a tip to a phylogenetic tree at a specified location.
#'
#' `AddTip()` extends \code{\link[ape]{bind.tree}}, which cannot handle
#'   single-taxon trees.
#'
#' @template treeParam
#' @param where The node or tip that should form the sister taxon to the new
#' node.  To add a new tip at the root, use `where = 0`.  By default, the
#' new tip is added to a random edge.
#' @param label Character string providing the label to apply to the new tip.
#' @param nodeLabel Character string providing a label to apply to the newly
#' created node, if `tree$node.label` is specified.
#' @param edgeLength Numeric specifying length of new edge. If `NULL`,
#' defaults to `lengthBelow`.
# Notice added in v1.10.0.9001, 2024-02-20:
#' This will become the default behaviour in a
#' future release; please manually specify the desired behaviour in your code.
#' @param lengthBelow Numeric specifying length below neighbour at which to
#' graft new edge. Values greater than the length of the edge will result
#' in negative edge lengths. If `NULL`, the default, the new tip will be added
#' at the midpoint of the broken edge. If inserting at the root (`where = 0`),
#' a new edge of length `lengthBelow` will be inserted.
#' If `NA`, the new leaf will be attached adjacent to `where`; at internal
#' nodes, this will result in polytomy.
#' 
#' @param nTip,nNode,rootNode Optional integer vectors specifying number of tips
#' and nodes in `tree`, and index of root node.
#' Not checked for correctness: specifying values here yields a marginal speed
#' increase at the cost of code safety.
#'
#' @return `AddTip()` returns a tree of class `phylo` with an additional tip
#' at the desired location.
#'
#' @template MRS
#'
#' @seealso Add one tree to another: \code{\link[ape]{bind.tree}()}
#'
#' @examples
#' tree <- BalancedTree(10)
#' 
#' # Add a leaf below an internal node
#' plot(tree)
#' ape::nodelabels() # Identify node numbers 
#' node <- 15        # Select location to add leaf
#' ape::nodelabels(bg = ifelse(NodeNumbers(tree) == node, "green", "grey"))
#'
#' plot(AddTip(tree, 15, "NEW_TIP"))
#' 
#' # Add edge lengths for an ultrametric tree
#' tree$edge.length <- rep(c(rep(1, 5), 2, 1, 2, 2), 2)
#' 
#' # Add a leaf to an external edge
#' leaf <- 5
#' plot(tree)
#' ape::tiplabels(bg = ifelse(seq_len(NTip(tree)) == leaf, "green", "grey"))
#' 
#' plot(AddTip(tree, 5, "NEW_TIP", edgeLength = NULL))
#' 
#' # Create a polytomy, rather than a new node
#' plot(AddTip(tree, 5, "NEW_TIP", edgeLength = NA))
#' 
#' @keywords tree
#' @family tree manipulation
#' @export
AddTip <- function(tree,
                   where = sample.int(tree[["Nnode"]] * 2 + 2L, size = 1) - 1L,
                   label = "New tip",
                   nodeLabel = "",
                   edgeLength = 0,
                   lengthBelow = NULL,
                   nTip = NTip(tree),
                   nNode = tree[["Nnode"]],
                   rootNode = RootNode(tree)
) {
  newTipNumber <- nTip + 1L
  treeEdge <- tree[["edge"]]
  edgeLengths <- tree[["edge.length"]]
  hasLengths <- !is.null(edgeLengths)
  
  if (is.character(where)) {
    tmp <- match(where, TipLabels(tree))
    if (is.na(tmp)) {
      stop("No tip labelled '", where, "'")
    }
    where <- tmp
  }
  
  # Find the row of "where" before renumbering
  if (where < 1L || where == rootNode) {
    case <- 1L
  } else {
    insertionEdge <- which(treeEdge[, 2] == where)
    case <- if (where <= nTip) 2L else 3L
  }
  # case = 1 -> y is bound on the root of x
  # case = 2 -> y is bound on a tip of x
  # case = 3 -> y is bound on a node of x
  
  # Because in all situations internal nodes need to be
  # renumbered, they are changed to negatives first, and
  # nodes eventually added will be numbered sequentially
  nodes <- treeEdge > nTip
  treeEdge[nodes] <- nTip - treeEdge[nodes]  # -1, ..., -nTip
  nextNode <- -nNode - 1L
  rootNode <- nTip - rootNode
  # Determine now, before we overwrite lengthBelow
  addingNode <- is.null(lengthBelow) || 
    (!is.null(lengthBelow) && !is.na(lengthBelow)) ||
    case == 2 
  
  switch(case, { # case = 1 -> y is bound on the root of x
    if (addingNode) {
      treeEdge <- rbind(c(nextNode, treeEdge[[1]]),
                        treeEdge,
                        c(nextNode, newTipNumber))
    } else {
      treeEdge <- rbind(treeEdge, c(rootNode, newTipNumber))
    }
    if (hasLengths) {
      if (is.null(lengthBelow) || is.na(lengthBelow)) {
        lengthBelow <- 0
      }
      edgeLengths <- c(if (addingNode) lengthBelow, edgeLengths,
                       if(is.null(edgeLength)) lengthBelow else edgeLength)
    }
    rootNode <- nextNode
  }, { # case = 2 -> y is bound on a tip of x
    beforeInsertion <- seq_len(insertionEdge)
    treeEdge[insertionEdge, 2] <- nextNode
    treeEdge <- rbind(treeEdge[beforeInsertion, ],
                      c(nextNode, where),
                      c(nextNode, newTipNumber),
                      treeEdge[-beforeInsertion, ])
    if (hasLengths) {
      if (is.null(lengthBelow)) {
        lengthBelow <- edgeLengths[insertionEdge] / 2L
      } else if (is.na(lengthBelow)) {
        lengthBelow <- 0
      }
      edgeLengths <- c(edgeLengths[beforeInsertion[-insertionEdge]],
                       edgeLengths[insertionEdge] - lengthBelow,
                       lengthBelow,
                       if(is.null(edgeLength)) lengthBelow else edgeLength,
                       edgeLengths[-beforeInsertion])
    }
  }, { # case = 3 -> y is bound on a node of x
    beforeInsertion <- seq_len(insertionEdge)
    
    if (addingNode) {
      treeEdge <- rbind(treeEdge[beforeInsertion, ],
                        c(nextNode, newTipNumber),
                        c(nextNode, treeEdge[insertionEdge, 2]),
                        treeEdge[-beforeInsertion, ])
      treeEdge[insertionEdge, 2] <- nextNode
    } else {
      treeEdge <- rbind(treeEdge[beforeInsertion, ],
                        c(treeEdge[insertionEdge, 2], newTipNumber),
                        treeEdge[-beforeInsertion, ])
    }
    
    if (hasLengths) {
      if (is.null(lengthBelow)) {
        lengthBelow <- edgeLengths[insertionEdge] / 2L
      } else if (is.na(lengthBelow)) {
        lengthBelow <- 0
      }
      edgeLengths <- if (addingNode) {
        c(edgeLengths[beforeInsertion[-insertionEdge]],
          edgeLengths[insertionEdge] - lengthBelow,
          if (is.null(edgeLength)) lengthBelow else edgeLength,
          lengthBelow, edgeLengths[-beforeInsertion])
      } else {
        c(edgeLengths[beforeInsertion], edgeLength,
          edgeLengths[-beforeInsertion])
      }
    }
    
  }
  )
  tree[["tip.label"]] <- c(tree[["tip.label"]], label)
  
  if (addingNode) {
    nNode <- nNode + 1L
    tree[["Nnode"]] <- nNode
  }
  
  ## renumber nodes:
  if (addingNode) {
    newNumbering <- integer(nNode)
    newNumbering[-rootNode] <- newTipNumber + 1L
    childNodes <- treeEdge[, 2] < 0L
    ## executed from right to left, so newNb is modified before x$edge:
    treeEdge[childNodes, 2] <-
      newNumbering[-treeEdge[childNodes, 2]] <-
      newTipNumber + 2:nNode
    treeEdge[, 1] <- newNumbering[-treeEdge[, 1]]
  } else {
    newNumbering <- newTipNumber + 1:nNode
    treeEdge[treeEdge < 0] <- newNumbering[-treeEdge[treeEdge < 0]]
  }
  
  
  
  tree[["edge"]] <- treeEdge
  if (hasLengths) {
    tree[["edge.length"]] <- edgeLengths
  }
  
  if (addingNode) {
    nodeLabels <- tree[["node.label"]]
    if (!is.null(nodeLabels)) {
      newLabels <- character(nNode)
      newLabels[newNumbering - newTipNumber] <- c(nodeLabels, "")
      tree[["node.label"]] <- newLabels
    }
  }
  
  # Return:
  tree
}

#' @rdname AddTip
#'
#' @details `AddTipEverywhere()` adds a tip to each edge in turn.
#'
#' @param includeRoot Logical; if `TRUE`, each position adjacent
#' to the root edge is considered to represent distinct edges; if `FALSE`,
#' they are treated as a single edge.
#' @return `AddTipEverywhere()` returns a list of class `multiPhylo` containing
#' the trees produced by adding `label` to each edge of `tree` in turn.
#'
#' @examples
#' # Set up multi-panel plot
#' oldPar <- par(mfrow = c(2, 4), mar = rep(0.3, 4), cex = 0.9)
#'
#' # Add leaf to each edge on a tree in turn
#' backbone <- BalancedTree(4)
#' # Treating the position of the root as instructive:
#' additions <- AddTipEverywhere(backbone, includeRoot = TRUE)
#' xx <- lapply(additions, plot)
#'
#' par(mfrow = c(2, 3))
#' # Don't treat root edges as distinct:
#' additions <- AddTipEverywhere(backbone, includeRoot = FALSE)
#' xx <- lapply(additions, plot)
#'
#' # Restore original plotting parameters
#' par(oldPar)
#'
#' @importFrom ape is.rooted
#' @export
AddTipEverywhere <- function(tree, label = "New tip", includeRoot = FALSE) {
  nTip <- NTip(tree)
  if (nTip == 0L) return(list(SingleTaxonTree(label)))
  if (nTip == 1L) return(list(StarTree(c(tree[["tip.label"]], label))))
  whichNodes <- seq_len(nTip + tree[["Nnode"]])
  edge <- tree[["edge"]]
  root <- RootNode(edge)
  if (!includeRoot) {
    parent <- edge[, 1]
    child <- edge[, 2]
    rootChildren <- child[parent == root]
    
    whichNodes <- if (length(rootChildren) == 2L && nTip > 2L) {
      rootChildrenNodes <- rootChildren[rootChildren > nTip]
      whichNodes[-c(root, rootChildrenNodes[1])]
    } else {
      whichNodes[-root]
    }
  }
  lapply(whichNodes, AddTip, tree = tree, label = label, nTip = nTip,
         rootNode = root)
}
