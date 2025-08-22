#' @importFrom methods new setClass
setClass("Splits", representation("raw"))

#' Convert object to `Splits`
#'
#' `as.Splits()` converts a phylogenetic tree to a `Splits` object representing
#' its constituent bipartition splits.
#'
#' @param x Object to convert into splits: perhaps a tree of class
#'  \code{\link[ape:read.tree]{phylo}}.
#'  If a logical matrix is provided, each row will be considered as a
#'  separate split.
#' @param tipLabels Character vector specifying sequence in which to order
#' tip labels.  Label order must (currently) match to combine or compare separate
#' `Splits` objects.
#' @param \dots Presently unused.
#'
#' @return `as.Splits()` returns an object of class `Splits`, or
#' (if `asSplits = FALSE`) a two-dimensional array of `raw` objects, with
#' each bit specifying whether or not the leaf corresponding to the respective
#' bit position is a member of the split.
#' Splits are named according to the node at the non-root end of the edge that
#' defines them. In rooted trees, the child of the rightmost root edge names
#' the split.
#'
#' @template MRS
#'
#' @examples
#' splits <- as.Splits(BalancedTree(letters[1:6]))
#' summary(splits)
#' TipsInSplits(splits)
#' summary(!splits)
#' TipsInSplits(!splits)
#'
#' length(splits + !splits)
#' length(unique(splits + !splits))
#'
#' summary(c(splits[[2:3]], !splits[[1:2]]))
#'
#' moreSplits <- as.Splits(PectinateTree(letters[6:1]), tipLabel = splits)
#' print(moreSplits, details = TRUE)
#' match(splits, moreSplits)
#' moreSplits %in% splits
#' 
#' as.Splits("....**", letters[1:6])
#'
#' @family Splits operations
#' @name Splits
#' @export
as.Splits <- function(x, tipLabels = NULL, ...) UseMethod("as.Splits")

#' @rdname Splits
#' @param asSplits Logical specifying whether to return a `Splits` object,
#'   or an unannotated two-dimensional array (useful where performance is
#'   paramount).
#' @export
as.Splits.phylo <- function(x, tipLabels = NULL, asSplits = TRUE, ...) {
  
  if (length(tipLabels) && !identical(x[["tip.label"]], tipLabels)) {
    x <- RenumberTips(x, tipLabels)
  }
  
  edge <- x[["edge"]]
  nEdge <- dim(edge)[[1]]
  order <- attr(x, "order")[[1]]
  edgeOrder <- if (length(order) == 0) {
    postorder_order(edge)
  } else {
    switch(order,
           "preorder" = rev(seq_len(nEdge)),
           "postorder" = seq_len(nEdge),
           postorder_order(edge))
  }

  # Return:
  edge_to_splits(edge, edgeOrder, tipLabels = x[["tip.label"]],
                 asSplits = asSplits, nTip = NTip(x), ...)
}

#' Efficiently convert edge matrix to splits
#' 
#' Wrapper for internal C++ function for maximum efficiency.
#' Improper input may crash R.  Behaviour not guaranteed.
#' It is advisable to contact the package maintainers before
#' relying on this function.
#' 
#' @template edgeParam
#' @param edgeOrder Integer vector such that `edge[edgeOrder, ]` returns a
#' postorder ordering of edges.
#' @param nTip Integer specifying number of leaves in tree.
#' @inheritParams Splits
#' @return `edge_to_splits()` uses the same return format as `as.Splits()`.
#' 
#' @seealso [`as.Splits()`][Splits] offers a safe access point to this
#' function that should be suitable for most users.
#' 
#' @family C++ wrappers
#' @export
edge_to_splits <- function(edge, edgeOrder, tipLabels = NULL, asSplits = TRUE,
                           nTip = NTip(edge), ...) {
  splits <- cpp_edge_to_splits(edge, edgeOrder - 1L, nTip)
  nSplits <- dim(splits)[[1]]

  # Return:
  if (asSplits) {
    structure(splits,
              nTip = nTip,
              tip.label = tipLabels,
              class = "Splits")
  } else {
    splits
  }
}

#' @rdname Splits
#' @export
as.Splits.multiPhylo <- function(x, tipLabels = unique(unlist(TipLabels(x))),
                                  asSplits = TRUE, ...) {
  lapply(x, as.Splits.phylo, tipLabels = tipLabels, asSplits = asSplits)
}

#' @rdname Splits
#' @export
as.Splits.Splits <- function(x, tipLabels = NULL, ...) {
  if (is.null(tipLabels)) {
    # Nothing needs doing.
    # Return:
    x
  } else {
    tipLabels <- TipLabels(tipLabels, single = TRUE)
    oldLabels <- attr(x, "tip.label")
    if (is.null(oldLabels)) {
      nTip <- attr(x, "nTip")
      if (length(tipLabels) == nTip) {
        attr(x, "tip.label") <- tipLabels
        x
      } else {
        stop(length(tipLabels), " labels provided; expecting ", nTip)
      }
    } else if (!identical(oldLabels, tipLabels)) {
      nTip <- attr(x, "nTip")
      if (length(x) == 0) {
        attr(x, "tip.label") <- tipLabels
        x
      } else {
        if (all(oldLabels %in% tipLabels) && all(tipLabels %in% oldLabels)) {
          as.Splits.logical(as.logical(x, tipLabels = tipLabels)
                            [, match(tipLabels, oldLabels), drop = FALSE],
                            tipLabels = tipLabels)
        } else {
          stop("Old and new labels must match")
        }
      }
    } else {
      x
    }
  }
}

#' @rdname Splits
#' @export
as.Splits.list <- function(x, tipLabels = NULL, asSplits = TRUE, ...) {
  if (inherits(x[[1]], "phylo")) {
    if (is.null(tipLabels)) {
      tipLabels <- unique(unlist(TipLabels(x)))
    }
    lapply(x, function(sp) as.Splits(
      sp,
      tipLabels = intersect(tipLabels, TipLabels(sp)),
      asSplits = asSplits
      ))
  } else if (inherits(x[[1]], "Splits")) {
    if (is.null(tipLabels)) {
      tipLabels <- attr(x, "tip.label")
      if (is.null(tipLabels)) {
        tipLabels <- attr(x[[1]], "tip.label")
      }
    }
    lapply(x, as.Splits, tipLabels = tipLabels, asSplits = asSplits)
  } else {
    stop("Unsupported list type; first item has class '",
         paste0(class(x[[1]]), collapse = ", "), "'")
  }
}

#' @rdname Splits
#' @export
as.Splits.matrix <- function(x, tipLabels = NULL, ...) {
  if (all(c("edge", "Nnode") %in% rownames(x))) {
    col1 <- x[, 1]
    if (is.list(col1)) {
      if (is.null(tipLabels)) {
        tipLabels <- col1[["tip.label"]]
        if (is.null(tipLabels)) {
          nTip <- dim(col1[["edge"]])[1] - col1[["Nnode"]] + 1L
          tipLabels <- seq_len(nTip)
          x <- rbind(x, replicate(ncol(x), list(tip.label = tipLabels)))
          rownames(x)[nrow(x)] <- "tip.label"
        }
      }
      lapply(seq_len(ncol(x)), function(i) {
        as.Splits(structure(x[, i], class = "phylo"), tipLabels = tipLabels)
      })
    } else {
      stop("Unsupported matrix. Columns should correspond to trees.")
    }
  } else if (is.numeric(x) && dim(x)[2] == 2) {
    edge_to_splits(x, postorder_order(x),
                   tipLabels = tipLabels, asSplits = TRUE, ...)
  } else {
    NextMethod()
  }
}

#' @rdname Splits
#' @export
as.Splits.logical <- function(x, tipLabels = NULL, ...) {
  dimX <- dim(x)
  
  if (is.null(dimX)) {
    nTip <- length(x)
    ret <- if (nTip == 0L) {
      matrix(raw(0), 0, 0)
    } else {
      pack_splits_logical_vec(x)
    }
  } else {
    nTip <- dimX[2L]
    ret <- `rownames<-`(pack_splits_logical(x), dimnames(x)[[1]])
  }
  
  if (is.null(tipLabels)) {
    # TODO when require R 4.1
    # tipLabels <- TipLabels(x) %||% TipLabels(nTip)
    tipLabels <- TipLabels(x)
    if (is.null(tipLabels)) {
      tipLabels <- TipLabels(nTip)
    }
  } else {
    tipLabels <- TipLabels(tipLabels)
  }
  
  structure(
    ret,
    nTip = nTip,
    tip.label = tipLabels,
    class = "Splits"
  )
}

#' @rdname Splits
#' @export
as.Splits.character <- function(x, tipLabels = NULL, ...) {
  nTip <- nchar(x[[1]])
  
  if (is.null(tipLabels)) {
    tipLabels <- TipLabels(x)
    if (length(tipLabels) == 1 && tipLabels == x) {
      tipLabels <- paste0("t", seq_len(nTip))
    }
  } else {
    tipLabels <- TipLabels(tipLabels)
  }
  
  sp <- .vapply(gregexpr("*", x, fixed = TRUE), tabulate, integer(nTip),
                nTip) > 0
  structure(t(.apply(sp, 2, as.Splits.logical)),
    nTip = nTip,
    tip.label = tipLabels,
    class = "Splits"
  )
}

#' @rdname Splits
#' @export
as.logical.Splits <- function(x, tipLabels = attr(x, "tip.label"), ...) {
  nTip <- attr(x, "nTip")
  if (dim(x)[1] == 0) {
    ret <- matrix(logical(0), 0, nTip)
  } else {
    ret <- matrix(as.logical(rawToBits(t(x))),
                  nrow = nrow(x), byrow = TRUE)[, seq_len(nTip), drop = FALSE]
  }
  dimnames(ret) <- list(rownames(x), tipLabels)
  ret
}

#' @family Splits operations
#' @export
print.Splits <- function(x, details = FALSE, ...) {
  nTip <- attr(x, "nTip")
  tipLabels <- attr(x, "tip.label")
  trivial <- TrivialSplits(x)
  cat(dim(x)[1], "bipartition", ifelse(dim(x)[1] == 1, "split", "splits"),
      if(any(trivial)) paste0("(", sum(trivial), " trivial)"),
      "dividing", nTip,
      if(is.null(tipLabels)) {
        "unlabelled tips."
      } else {
        if (nTip) {
          if (nTip == 1) {
            paste("tip,", tipLabels[1])
          } else {
            paste("tips,", tipLabels[1], "..", tipLabels[nTip])
          }
        } else {
          "tips"
        }
      }
      )
  if (details) {
    splitNames <- rownames(x)
    if (!is.null(splitNames)) {

      nameLengths <- vapply(splitNames, nchar, 0)
      namePads <- vapply(nameLengths, function(thisLength)
        paste0(rep.int(" ", max(nameLengths) - thisLength), collapse=""), character(1))
      splitNames <- paste0(splitNames, namePads)
    } else {
      splitNames <- character(length(x))
      nameLengths = 0L
    }
    cat("\n ", paste0(rep.int(" ", max(nameLengths)), collapse = ""),
        paste0(rep_len(c(1:9, " "), nTip), collapse = ""))

    for (i in seq_len(dim(x)[1])) {
      split <- x[i, , drop = FALSE]
      cat("\n", splitNames[i], "",
          paste(ifelse(as.logical(rawToBits(split)[seq_len(nTip)]), "*", "."),
                collapse = ""))
    }
  }
}

#' @family Splits operations
#' @importFrom utils head
#' @export
head.Splits <- function(x, n = 6L, ...) {
  if (n < 0) {
    n <- max(0L, length(x) + n)
  }
  if (length(x) > n) {
    x[[seq_len(n)]]
  } else {
    x
  }
}

#' @family Splits operations
#' @importFrom utils tail
#' @export
tail.Splits <- function(x, n = 6L, ...) {
  if (n < 0) {
    n <- max(0L, length(x) + n)
  }
  if (length(x) > n) {
    x[[length(x) + 1L - rev(seq_len(n))]]
  } else {
    x
  }
}

#' @family Splits operations
#' @export
summary.Splits <- function(object, ...) {
  print(object, details = TRUE, ...)
  nTip <- attr(object, "nTip")
  if (is.null(attr(object, "tip.label"))) {
    cat("\n\nTips not labelled.")
  } else {
    cat("\n\n", paste0("Tip ", seq_len(nTip), ": ", attr(object, "tip.label"),
                     "\t", c(character(4L), "\n")[seq_len(min(nTip, 5L))]))
  }
}

#' @family Splits operations
#' @export
names.Splits <- function(x) rownames(x)


#' @family Splits operations
#' @export
`names<-.Splits` <- `rownames<-`


#' @family Splits operations
#' @importFrom stringi stri_paste
#' @export
as.character.Splits <- function(x, ...) {
  tipLabels <- attr(x, "tip.label")
  nTip <- attr(x, "nTip")

  apply(as.logical(x), 1L, function(inSplit) {
    stri_paste(stri_paste(tipLabels[inSplit], collapse = " "), " | ",
               stri_paste(tipLabels[!inSplit], collapse = " "))
  })

}

#' @family Splits operations
#' @export
as.phylo.Splits <- function(x, ...) {
  .PreorderTree(
    edge = splits_to_edge(x, NTip(x)),
    tip.label = TipLabels(x)
  )
}

#' @family Splits operations
#' @export
c.Splits <- function(...) {
  splits <- list(...)
  nTip <- unique(vapply(splits, attr, 1, "nTip"))
  if (length(nTip) > 1L) {
    stop("Splits must relate to identical tips.")
  }
  tips <- lapply(splits, attr, "tip.label")
  if (length(unique(lapply(tips, sort))) > 1L) {
    stop("All splits must bear identical tips")
  }
  tipLabels <- tips[[1]]
  splits <- c(splits[1], lapply(splits[seq_along(splits)[-1]], as.Splits,
                                  tipLabels = tipLabels))

  structure(do.call(rbind, splits),
            nTip = nTip,
            tip.label = tipLabels,
            class = "Splits")
}

#' @family Splits operations
#' @export
`+.Splits` <- function(...) c(...)


#' @family Splits operations
#' @export
`-.Splits` <- function(x, y) {
  x[[!x %in% y]]
}

#' @family Splits operations
#' @export
`[[.Splits` <- function(x, ..., drop = TRUE) {
  elements <- list(...)
  if (length(elements) == 0L) {
    x[]
  } else if (length(elements) == 1L) {
    ret <- x[elements[[1]], , drop = FALSE]
    at <- attributes(x)
    at[["dim"]] <- dim(ret)
    at[["dimnames"]][[1]] <- rownames(x)[elements[[1]]]
    attributes(ret) <- at
    ret
  } else {
    stop("Too many dimensions specified")
  }
}

.MaskSplits <- function(x) {
  mask_splits(x)
}

#' @family Splits operations
#' @method ! Splits
#' @export
`!.Splits` <- function(x) {
  not_splits(x)
}

#' @family Splits operations
#' @method & Splits
#' @export
`&.Splits` <- function(e1, e2) {
  and_splits(e1, e2)
}

#' @family Splits operations
#' @method | Splits
#' @export
`|.Splits` <- function(e1, e2) {
  or_splits(e1, e2)
}


#' Exclusive OR operation
#' @rdname xor
#' @param x,y Objects to be compared.
setGeneric("xor")

#' @family Splits operations
#' @rdname xor
#' @importFrom methods setMethod
#' @export
setMethod("xor", signature = representation(x = "Splits", y = "Splits"),
          function(x, y) xor_splits(x, y))

#' @export
t.Splits <- function(x) t(x[])

#' @family Splits operations
#' @export
length.Splits <- function(x) nrow(x)

#' @family Splits operations
#' @importFrom stats setNames
#' @export
duplicated.Splits <- function(x, incomparables, fromLast = FALSE,
                              withNames = TRUE, ...) {
  if (missing(incomparables)) {
    ret <- duplicated_splits(x, isTRUE(fromLast))
    if (withNames) names(ret) <- names(x)
    ret
  } else {
    dupX <- !x
    useOrig <- x[, 1] < dupX[, 1]
    dupX[useOrig, ] <- x[useOrig, ]
    
    if (dim(dupX)[2] == 1) {
      ret <- duplicated.default(dupX, incomparables = incomparables, 
                         fromLast = fromLast, ...)
      if (withNames) {
        setNames(ret, names(dupX))
      } else {
        ret
      }
    } else {
      ret <- duplicated.array(dupX, MARGIN = 1, incomparables = incomparables,
                              fromLast = fromLast, ...)
      if (withNames) {
        ret
      } else {
        as.logical(ret)
      }
    }
  }
}

#' @family Splits operations
#' @export
unique.Splits <- function(x, incomparables, ...) {
  at <- attributes(x)
  dups <- duplicated(x, incomparables, ...)
  x <- x[!dups, , drop = FALSE]
  at[["dim"]] <- dim(x)
  at[["dimnames"]] <- dimnames(x)
  attributes(x) <- at
  x
}

#' @family Splits operations
#' @export
rev.Splits <- function(x) {
  newOrder <- seq.int(dim(x)[1], 1)
  x[] <- x[newOrder, ]
  rownames(x) <- rownames(x)[newOrder]
  x
}

#' Polarize splits on a single taxon
#'
#' @param x Object of class [`Splits`].
#' @param pole Numeric or character identifying tip that should polarize each
#' split.
#'
#' @return `PolarizeSplits()` returns a `Splits` object in which `pole` is
#' represented by a zero bit
#' @family Splits operations
#' @export
PolarizeSplits <- function(x, pole = 1L) {
  nTip <- attr(x, "nTip")
  if (!is.numeric(pole)) {
    pole <- match(pole, attr(x, "tip.label"))
  }
  if (is.na(pole)) {
    stop("Could not find `pole` in labels of `x`")
  }
  if (pole < 1) {
    stop("`pole` must be positive")
  }
  if (pole > nTip) {
    stop("`pole` must correspond to a leaf")
  }

  poleMask <- logical(8 * ceiling(nTip / 8))
  poleMask[pole] <- T
  poleMask <- packBits(poleMask)
  flip <- !apply(t(x) & poleMask, 2, function(x) any(as.logical(x)))

  x[flip, ] <- !x[flip, ]
  remainder <- (8L - nTip) %% 8L
  if (remainder) {
    lastSplit <- dim(x)[2]
    endMask <- packBits(c(!logical(8L - remainder), logical(remainder)))
    x[, lastSplit] <- x[, lastSplit] & endMask
  }
  x
}
