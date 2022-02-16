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
#' dropped, or integer vector specifying the indices of leaves to be dropped.
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
#' plot(DropTip(tree, c('t5', 't6')))
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
#' @export
DropTip.phylo <- function(tree, tip, preorder = TRUE, check = TRUE) {
  if (preorder) {
    tree <- Preorder(tree)
  }
  labels <- tree[["tip.label"]]
  nTip <- length(labels)
  if (is.null(tip) || !length(tip) || any(is.na(tip))) {
    drop <- logical(0)
  } else if (is.character(tip)) {
    drop <- match(tip, labels)
    if (check) {
      missing <- is.na(drop)
      if (any(missing)) {
        warning(paste(tip[missing], collapse = ', '), " not present in tree")
        drop <- drop[!missing]
      }
    }
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
        drop <- c(tip[tip <= nTip],
                  .DescendantTips(parent, child, nTip, tip[tip > nTip]))
      } else {
        drop <- tip
      }
    } else {
      drop <- tip
    }
  } else {
    stop("`tip` must be of type character or numeric")
  }
  
  nDrop <- length(drop)
  if (nDrop) {
    if (nDrop == nTip) {
      return(structure(list(edge = matrix(0, 0, 2), tip.label = character(0),
                            Nnode = 0), class = 'phylo'))
    }
    
    keep <- !tabulate(drop, nbins = nTip)
    weights <- tree[["edge.length"]]
    if (!is.null(weights)) {
      original <- PathLengths(tree)
      rownames(original) <- apply(original[, 1:2], 1, paste, collapse = " ")
      verts <- KeptVerts(tree, keep)
    }
    tree[["edge"]] <- keep_tip(tree[["edge"]], keep)
    tree[["tip.label"]] <- labels[-drop]
    tree[["Nnode"]] <- dim(tree[["edge"]])[1] + 1L - (nTip - nDrop)
    
    if (!is.null(weights)) {
      tree[["edge.length"]] <- original[
        apply(matrix(which(verts)[tree[["edge"]]], ncol = 2L),
              1, paste, collapse = " "), "length"]
    }
  }
  
  # Return:
  tree
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
    attr(tree, 'TipLabel') <- setdiff(at[["TipLabel"]], tip)
  }
  tree
}

#' @rdname DropTip
#' @return `KeepTip()` returns `tree` with all leaves not in `tip` removed,
#' in preorder.
#' @export
KeepTip <- function(tree, tip, preorder = TRUE, check = TRUE) {
  labels <- if (is.character(tip)) {
    TipLabels(tree)
  } else {
    seq_len(NTip(tree))
  }
  if (!all(tip %in% labels)) {
    warning("Tips not in tree: ", paste0(setdiff(tip, labels), collapse = ', '))
  }
  keep <- setdiff(labels, tip)
  DropTip(tree, keep, preorder, check)
}
