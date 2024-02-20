#' Add a tip to a phylogenetic tree
#'
#' `AddTip()` adds a tip to a phylogenetic tree at a specified location.
#'
#' `AddTip()` extends \code{\link{bind.tree}}, which cannot handle
#'   single-taxon trees.
#'
#' @template treeParam
#' @param where The node or tip that should form the sister taxon to the new
#' node.  To add a new tip at the root, use `where = 0`.  By default, the
#' new tip is added to a random edge.
#' @param label Character string providing the label to apply to the new tip.
#' @param edgeLength Numeric specifying length of new edge
#' @param lengthBelow Numeric specifying length below neighbour at which to
#' graft new edge. Values greater than the length of the edge will result
#' in negative edge lengths. If `NULL`, the default, the new tip will be added
#' at the midpoint of the broken edge. If inserting at the root (`where = 0`),
#' a new edge of length `lengthBelow` will be inserted.
#' @param nTip,nNode,rootNode Optional integer vectors specifying number of tips and
#' nodes in `tree`, and index of root node.
#' Not checked for correctness: specifying values here trades code safety for a
#' nominal speed increase.
#'
#' @return `AddTip()` returns a tree of class `phylo` with an additional tip
#' at the desired location.
#'
#' @template MRS
#'
#' @seealso Add one tree to another: \code{\link{bind.tree}()}
#'
#' @examples
#' tree <- BalancedTree(10)
#' 
#' # Add a leaf below an internal node
#' plot(tree)
#' ape::nodelabels()
#' ape::nodelabels(15, 15, bg = "green")
#'
#' plot(AddTip(tree, 15, "NEW_TIP"))
#' 
#' # Add a leaf to an external edge
#' plot(tree)
#' ape::tiplabels()
#' ape::tiplabels(5, 5, bg = "green")
#' 
#' plot(AddTip(tree, 5, "NEW_TIP"))
#' 
#' @keywords tree
#' @family tree manipulation
#'
#' @export
AddTip <- function(tree,
                   where = sample.int(tree[["Nnode"]] * 2 + 2L, size = 1) - 1L,
                   label = "New tip",
                   edgeLength = 0,
                   lengthBelow = NULL,
                   nTip = NTip(tree),
                   nNode = tree[["Nnode"]],
                   rootNode = RootNode(tree)
) {
  newTipNumber <- nTip + 1L
  treeEdge <- tree[["edge"]]
  edgeLengths <- tree[["edge.length"]]
  lengths <- !is.null(edgeLengths)
  
  if (is.character(where)) {
    tmp <- match(where, TipLabels(tree))
    if (is.na(tmp)) stop("No tip labelled '", where, "'")
    where <- tmp
  }
  ## find the row of "where" before renumbering
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
  
  switch(case, { # case = 1 -> y is bound on the root of x
    treeEdge <- rbind(c(nextNode, treeEdge[1]), treeEdge, c(nextNode, newTipNumber))
    if (lengths) {
      if (is.null(lengthBelow)) {
        lengthBelow <- 0
      }
      edgeLengths <- c(lengthBelow, edgeLengths, edgeLength)
    }
    rootNode <- nextNode
  }, { # case = 2 -> y is bound on a tip of x
    beforeInsertion <- seq_len(insertionEdge)
    treeEdge[insertionEdge, 2] <- nextNode
    treeEdge <- rbind(treeEdge[beforeInsertion, ],
                      c(nextNode, where),
                      c(nextNode, newTipNumber),
                      treeEdge[-beforeInsertion, ])
    if (lengths) {
      if (is.null(lengthBelow)) {
        lengthBelow <- edgeLengths[insertionEdge] / 2L
      }
      edgeLengths <- c(edgeLengths[beforeInsertion[-insertionEdge]],
                       edgeLengths[insertionEdge] - lengthBelow,
                       lengthBelow,
                       edgeLength,
                       edgeLengths[-beforeInsertion])
    }
  }, { # case = 3 -> y is bound on a node of x
    beforeInsertion <- seq_len(insertionEdge)
    
    treeEdge <- rbind(treeEdge[beforeInsertion, ],
                      c(nextNode, newTipNumber),
                      c(nextNode, treeEdge[insertionEdge, 2]),
                      treeEdge[-beforeInsertion, ])
    treeEdge[insertionEdge, 2] <- nextNode
    
    if (lengths) {
      if (is.null(lengthBelow)) {
        lengthBelow <- edgeLengths[insertionEdge] / 2L
      }
      edgeLengths <- c(edgeLengths[beforeInsertion[-insertionEdge]],
                       edgeLengths[insertionEdge] - lengthBelow,
                       edgeLength,
                       lengthBelow,
                       edgeLengths[-beforeInsertion])
    }
    
  }
  )
  tree[["tip.label"]] <- c(tree[["tip.label"]], label)
  
  nNode <- nNode + 1L
  tree[["Nnode"]] <- nNode
  
  ## renumber nodes:
  newNumbering <- integer(nNode)
  newNumbering[-rootNode] <- newTipNumber + 1L
  childNodes <- treeEdge[, 2] < 0L
  
  ## executed from right to left, so newNb is modified before x$edge:
  treeEdge[childNodes, 2] <-
    newNumbering[-treeEdge[childNodes, 2]] <-
    newTipNumber + 2:nNode
  treeEdge[, 1] <- newNumbering[-treeEdge[, 1]]
  
  tree[["edge"]] <- treeEdge
  if (lengths) {
    tree[["edge.length"]] <- edgeLengths
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
#' oldPar <- par(mfrow = c(2, 4), mar = rep(0.3, 4), cex = 0.9)
#'
#' backbone <- BalancedTree(4)
#' # Treating the position of the root as instructive:
#' additions <- AddTipEverywhere(backbone, includeRoot = TRUE)
#' xx <- lapply(additions, plot)
#'
#' par(mfrow=c(2, 3))
#' # Don't treat root edges as distinct:
#' additions <- AddTipEverywhere(backbone, includeRoot = FALSE)
#' xx <- lapply(additions, plot)
#'
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
