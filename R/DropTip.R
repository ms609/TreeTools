#' Drop leaves from tree
#'
#' `DropTip()` removes specified leaves from a phylogenetic tree, collapsing
#' incident branches.
#'
#' This function differs from [`ape::drop.tip()`], which roots unrooted trees,
#' and which can crash when trees' internal numbering follows unexpected schema.
#'
#' @template treeParam
#' @param tip Character vector specifying labels of leaves in tree to be
#' dropped or kept, or integer vector specifying the indices of leaves to be
#' dropped or kept.
#' Specifying the index of an internal node will drop all descendants of that
#' node.
#' @param preorder Logical specifying whether to [Preorder] `tree` before
#' dropping tips.  Specifying `FALSE` saves a little time, but will result in
#' undefined behaviour if `tree` is not in preorder.
#' @param check Logical specifying whether to check validity of `tip`. If
#' `FALSE` and `tip` contains entries that do not correspond to leaves of the
#' tree, undefined behaviour may occur.
#'
#' @return `DropTip()` returns a tree of class `phylo`, with the requested
#' leaves removed. The edges of the tree will be numbered in preorder,
#' but their sequence may not conform to the conventions of [`Preorder()`].
#'
#' @examples
#' tree <- BalancedTree(9)
#' plot(tree)
#' plot(DropTip(tree, c("t5", "t6")))
#' 
#' unrooted <- UnrootTree(tree)
#' plot(unrooted)
#' plot(DropTip(unrooted, 4:5))
#'
#' @family tree manipulation
#' @template MRS
#' @export
DropTip <- function(tree, tip, preorder = TRUE, check = TRUE) {
  UseMethod("DropTip")
}

#' @rdname DropTip
#' @return `KeepTip()` returns `tree` with all leaves not in `tip` removed,
#' in preorder.
#' @export
KeepTip <- function(tree, tip, preorder = TRUE, check = TRUE) {
  UseMethod("KeepTip")
}

#' @rdname DropTip
#' @export
DropTip.phylo <- function(tree, tip, preorder = TRUE, check = TRUE) {
  if (preorder && (length(attr(tree, "order")) == 0 || 
                   attr(tree, "order")[[1]] != "preorder")) {
    tree <- Preorder(tree)
  }
  labels <- tree[["tip.label"]]
  nTip <- length(labels)
  
  # Establish which tips should be dropped
  if (is.null(tip) || !length(tip) || any(is.na(tip))) {
    drop <- logical(nTip)
  } else if (is.character(tip)) {
    if (check) {
      if (sum(unique(labels) %in% tip) != length(tip)) {
        warning(paste(tip[!tip %in% labels], collapse = ", "),
                " not present in tree")
      }
    }
    drop <- labels %in% tip
  } else if (is.numeric(tip)) {
    if (check) {
      nNodes <- nTip + tree[["Nnode"]]
      if (any(tip > nNodes)) {
        warning("Tree only has ", nNodes, " nodes")
        tip <- tip[tip <= nNodes]
      }
      if (any(tip < 1L)) {
        warning("`tip` must be > 0")
        tip <- tip[tip > 0L]
      }
      
      if (any(tip > nTip)) {
        edge <- tree[["edge"]]
        parent <- edge[, 1]
        child <- edge[, 2]
        drop <- tabulate(c(tip[tip <= nTip],
                           .DescendantTips(parent, child, nTip,
                                           tip[tip > nTip])),
                         nbins = nTip)
      } else {
        drop <- as.logical(tabulate(tip, nbins = nTip))
      }
    } else {
      drop <- as.logical(tabulate(tip, nbins = nTip))
    }
  } else if (is.logical(tip)) {
    if (check) {
      if (length(tip) != nTip) {
        stop("`tip` must list `TRUE` or `FALSE` for each leaf.")
      }
    }
    drop <- tip
  } else {
    stop("`tip` must be of type character or numeric")
  }
  
  # Having established which tips to drop,
  # work out what to do about it
  if (all(drop)) {
    return(structure(list(edge = matrix(0, 0, 2), tip.label = character(0),
                          Nnode = 0), class = "phylo"))
  }
  
  if (any(drop)) {
    weights <- tree[["edge.length"]]
    keep <- !drop
    nodeLabel <- tree[["node.label"]]
    if (!is.null(weights)) {
      original <- PathLengths(tree, fullMatrix = TRUE)
    }
    if (!is.null(weights) || !is.null(nodeLabel)) {
      verts <- KeptVerts(tree, keep)
      tree[["node.label"]] <- nodeLabel[verts[-seq_len(nTip)]]
    }
    tree[["edge"]] <- keep_tip(tree[["edge"]], keep)
    tree[["tip.label"]] <- labels[keep]
    tree[["Nnode"]] <- dim(tree[["edge"]])[1] + 1L - sum(keep)
    
    if (!is.null(weights)) {
      newVerts <- matrix(which(verts)[tree[["edge"]]], ncol = 2L)
      tree[["edge.length"]] <- original[newVerts]
    }
  }
  
  # Return:
  tree
}

#' @rdname DropTip
#' @export
KeepTip.phylo <- function(tree, tip, preorder = TRUE, check = TRUE) {
  drop <- switch(
    mode(tip),
    "logical" = {
      !tip
    }, "character" = {
      setdiff(TipLabels(tree), MatchStrings(tip, TipLabels(tree), warning))
    }, "numeric" = {
      nTip <- NTip(tree)
      nValid <- if (inherits(tree, "phylo")) {
        nTip + tree[["Nnode"]]
      } else {
        nTip
      }
      badTips <- tip > nValid | tip < 1
      if (any(badTips)) {
        warning("`tip` must be between 1 and ", nValid)
        tip <- tip[!badTips]
      }
      nodes <- tip > nTip
      keep <- if (any(nodes)) {
        edge <- tree[["edge"]]
        theNodes <- tip[nodes]
        nodeDesc <- DescendantTips(edge[, 1], edge[, 2], node = theNodes)
        union(tip[!nodes], which(
          if (is.null(dim(nodeDesc))) nodeDesc else colSums(nodeDesc) > 0))
      } else {
        tip
      }
      setdiff(seq_len(nTip), keep)
    },
    stop("`tip` must be of type character, logical or numeric")
  )
  DropTip(tree, drop, preorder, check)
}

#' @rdname DropTip
#' @export
KeepTip.Splits <- KeepTip.phylo

#' @rdname DropTip
#' @examples summary(DropTip(as.Splits(tree), 4:5))
#' @family split manipulation functions
#' @export
DropTip.Splits <- function(tree, tip, preorder, check = TRUE) {
  labels <- TipLabels(tree)
  if (is.null(tip) || !length(tip) || any(is.na(tip))) {
    return(tree)
  } else if (is.character(tip)) {
    if (check) {
      drop <- match(tip, labels)
      missing <- is.na(drop)
      if (any(missing)) {
        warning(paste(tip[missing], collapse = ", "), " not present in tree")
        drop <- drop[!missing]
      }
      drop <- as.logical(tabulate(drop, NTip(tree)))
    } else {
      drop <- labels %in% tip
    }
  } else if (is.numeric(tip)) {
    nTip <- NTip(tree)
    if (check) {
      if (any(tip > nTip)) {
        warning("Tree only has ", nTip, " leaves")
        tip <- tip[tip <= nTip]
      }
      if (any(tip < 1L)) {
        warning("`tip` must be > 0")
        tip <- tip[tip > 0L]
      }
    }
    drop <- as.logical(tabulate(tip, nTip))
  } else if (is.logical(tip)) {
    if (check && length(tip) != NTip(tree)) {
      stop("`tip` must contain an entry for each leaf of `tree`.")
    }
    drop <- tip
  } else {
    stop("`tip` must be of type character, logical or numeric")
  }
  
  keep <- !drop
  thinner <- structure(thin_splits(tree, drop),
                       nTip = sum(keep),
                       tip.label = labels[keep],
                       class = "Splits")
  
  # Return:
  unique(thinner, fromLast = TRUE)
}

# nodes must all be internal
.DescendantTips <- function(parent, child, nTip, nodes, isDesc = logical(nTip)) {
  newDescs <- child[parent %in% nodes]
  recurse <- newDescs > nTip
  
  # Return:
  if (any(recurse)) {
    isDesc[newDescs[!recurse]] <- TRUE
    .DescendantTips(parent, child, nTip, newDescs[recurse], isDesc)
  } else {
    isDesc[newDescs] <- TRUE
    which(isDesc)
  }
}

#' @describeIn DropTip Direct call to `DropTip.phylo()`, to avoid overhead of
#' querying object's class.
#' @export
DropTipPhylo <- DropTip.phylo

#' @rdname DropTip
#' @export
DropTip.multiPhylo <- function(tree, tip, preorder = TRUE, check = TRUE) {
  at <- attributes(tree)
  tree <- lapply(tree, DropTip, tip, preorder)
  attributes(tree) <- at
  if (!is.null(at[["TipLabel"]])) {
    attr(tree, "TipLabel") <- setdiff(at[["TipLabel"]], tip)
  }
  tree
}

#' @rdname DropTip
#' @export
KeepTip.multiPhylo <- function(tree, tip, preorder = TRUE, check = TRUE) {
  at <- attributes(tree)
  tree <- lapply(tree, KeepTip, tip = tip, preorder = preorder, check = check)
  attributes(tree) <- at
  if (!is.null(at[["TipLabel"]])) {
    attr(tree, "TipLabel") <- intersect(at[["TipLabel"]], tip)
  }
  tree
}

#' @rdname DropTip
#' @export
DropTip.list <- function(tree, tip, preorder = TRUE, check = TRUE) {
  if (all(vapply(tree, inherits, logical(1), "phylo"))) {
    DropTip.multiPhylo(tree, tip, preorder, check)
  } else {
    NextMethod()
  }
}

#' @rdname DropTip
#' @export
KeepTip.list <- function(tree, tip, preorder = TRUE, check = TRUE) {
  if (all(vapply(tree, inherits, logical(1), "phylo"))) {
    KeepTip.multiPhylo(tree, tip, preorder, check)
  } else {
    NextMethod()
  }
}

#' @rdname DropTip
#' @export
DropTip.NULL <-  function(tree, tip, preorder = TRUE, check = TRUE) {
  NULL
}

#' @rdname DropTip
#' @export
KeepTip.NULL <-  function(tree, tip, preorder = TRUE, check = TRUE) {
  NULL
}

#' @describeIn DropTip Faster version with no checks.
#' Does not retain labels or edge weights.
#' Edges must be listed in preorder.
#' May crash if improper input is specified.
#' @export
KeepTipPreorder <- function(tree, tip) {
  if (!any(tip)) {
    structure(list(edge = matrix(0, 0, 2), tip.label = character(0),
                          Nnode = 0), class = "phylo")
  } else if (all(tip)) {
    tree
  } else {
    tree[["edge"]] <- keep_tip(tree[["edge"]], tip)
    tree[["tip.label"]] <- as.character(which(tip))
    tree[["Nnode"]] <- dim(tree[["edge"]])[1] + 1L - sum(tip) 
    tree
  }
}
  
#' @describeIn DropTip Faster version with no checks.
#' Does not retain labels or edge weights.
#' Edges must be listed in postorder.
#' May crash if improper input is specified.
#' @export
KeepTipPostorder <- function(tree, tip) {
  if (!any(tip)) {
    structure(list(edge = matrix(0, 0, 2), tip.label = character(0),
                          Nnode = 0), class = "phylo")
  } else if (all(tip)) {
    tree
  } else {
    edge <- tree[["edge"]]
    tree[["edge"]] <- keep_tip(edge[dim(edge)[1]:1, ], tip)
    tree[["tip.label"]] <- tree[["tip.label"]][tip]
    tree[["Nnode"]] <- dim(tree[["edge"]])[1] + 1L - sum(tip) 
    tree
  }
}
