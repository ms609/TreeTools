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
#' @family Splits operations
#' @name Splits
#' @export
as.Splits <- function(x, tipLabels = NULL, ...) UseMethod('as.Splits')

#' @rdname Splits
#' @param asSplits Logical specifying whether to return a `Splits` object,
#'   or an unannotated two-dimensional array (useful where performance is
#'   paramount).
#' @export
as.Splits.phylo <- function(x, tipLabels = NULL, asSplits = TRUE, ...) {
  if (!is.null(tipLabels)) {
    x <- RenumberTips(x, tipLabels)
  }
  edge <- x[["edge"]]
  nEdge <- dim(edge)[1]
  order <- attr(x, "order")[1]
  edgeOrder <- if (is.null(order)) {
    postorder_order(edge)
  } else {
    switch(order,
           "preorder" = nEdge:1,
           "postorder" = seq_len(nEdge),
           postorder_order(edge))
  }

  # Return:
  .as.Splits.edge(edge, edgeOrder, tipLabels = x[["tip.label"]],
                  asSplits = asSplits, nTip = NTip(x), ...)
}

.as.Splits.edge <- function(edge, edgeOrder, tipLabels = NULL, asSplits = TRUE,
                             nTip = NTip(edge), ...) {
  splits <- cpp_edge_to_splits(edge, edgeOrder - 1L, nTip)
  nSplits <- dim(splits)[1]

  # Return:
  if (asSplits) {
    structure(splits,
              nTip = nTip,
              tip.label = tipLabels,
              class = 'Splits')
  }
  else {
    splits
  }
}

#' @rdname Splits
#' @export
as.Splits.multiPhylo <- function(x, tipLabels = x[[1]][["tip.label"]],
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
    oldLabels <- attr(x, 'tip.label')
    if (is.null(oldLabels)) {
      nTip <- attr(x, 'nTip')
      if (length(tipLabels) == nTip) {
        attr(x, 'tip.label') <- tipLabels
        x
      } else {
        stop(length(tipLabels), " labels provided; expecting ", nTip)
      }
    } else if (!identical(oldLabels, tipLabels)) {
      nTip <- attr(x, 'nTip')
      if (length(x) == 0) {
        attr(x, 'tip.label') <- tipLabels
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
  if (inherits(x[[1]], 'phylo')) {
    if (is.null(tipLabels)) {
      tipLabels <- x[[1]][["tip.label"]]
    }
    lapply(x, as.Splits, tipLabels = tipLabels, asSplits = asSplits)
  } else if (inherits(x[[1]], 'Splits')) {
    if (is.null(tipLabels)) {
      tipLabels <- attr(x, 'tip.label')
      if (is.null(tipLabels)) {
        tipLabels <- attr(x[[1]], 'tip.label')
      }
    }
    lapply(x, as.Splits, tipLabels = tipLabels, asSplits = asSplits)
  } else {
    stop("Unsupported list type.")
  }
}

#' @rdname Splits
#' @export
as.Splits.matrix <- function(x, tipLabels = NULL, ...) {
  if (all(c('edge', 'Nnode') %in% rownames(x))) {
    col1 <- x[, 1]
    if (is.list(col1)) {
      if (is.null(tipLabels)) {
        tipLabels <- col1[["tip.label"]]
        if (is.null(tipLabels)) {
          nTip <- dim(col1[["edge"]])[1] - col1[["Nnode"]] + 1L
          tipLabels <- seq_len(nTip)
          x <- rbind(x, replicate(ncol(x), list(tip.label = tipLabels)))
          rownames(x)[nrow(x)] <- 'tip.label'
        }
      }
      lapply(seq_len(ncol(x)), function(i) {
        as.Splits(structure(x[, i], class = 'phylo'), tipLabels = tipLabels)
      })
    } else {
      stop("Unsupported matrix. Columns should correspond to trees.")
    }
  } else if (dim(x)[2] == 2) {
    .as.Splits.edge(x, postorder_order(x),
                    tipLabels = NULL, asSplits = TRUE, ...)
  } else {
    NextMethod()
  }
}

#' @rdname Splits
#' @export
as.Splits.logical <- function(x, tipLabels = NULL, ...) {
  powersOf2 <- as.raw(c(1L, 2L, 4L, 8L, 16L, 32L, 64L, 128L))
  dimX <- dim(x)
  if (is.null(dimX)) {
    nTip <- length(x)

    if (is.null(tipLabels)) {
      tipLabels <- TipLabels(x)
      if (is.null(tipLabels)) {
        tipLabels <- paste0('t', seq_len(nTip))
      }
    } else {
      tipLabels <- TipLabels(tipLabels)
    }

    structure(matrix(packBits(c(x, logical((8L - nTip) %% 8))), nrow = 1L),
              nTip = nTip,
              tip.label = tipLabels,
              class = 'Splits')
  } else {
    if (is.null(tipLabels)) {
      tipLabels <- TipLabels(x)
    }
    nTip <- dimX[2]
    if (is.null(tipLabels)) {
      tipLabels <- paste0('t', seq_len(nTip))
    }

    nBin <- (nTip %% 8 != 0) + (nTip / 8)
    structure(
      matrix(packBits(t(cbind(x, matrix(F, dimX[1], (8L - nTip) %% 8)))),
             nrow = dimX[1], ncol = nBin,
             byrow = TRUE, dimnames = list(rownames(x), NULL)),
      nTip = nTip,
      tip.label = tipLabels,
      class = 'Splits'
    )
  }
}

#' @rdname Splits
#' @export
as.logical.Splits <- function(x, tipLabels = NULL, ...) {
  nTip <- attr(x, 'nTip')
  if (dim(x)[1] == 0) {
    ret <- matrix(logical(0), 0, nTip)
  } else {
    ret <- matrix(as.logical(rawToBits(t(x))),
                  nrow = nrow(x), byrow = TRUE)[, seq_len(nTip), drop = FALSE]
  }
  dimnames(ret) <- list(rownames(x), attr(x, 'tip.label'))
  ret
}

#' @family Splits operations
#' @export
print.Splits <- function(x, details = FALSE, ...) {
  nTip <- attr(x, 'nTip')
  tipLabels <- attr(x, 'tip.label')
  trivial <- TrivialSplits(x)
  cat(dim(x)[1], "bipartition", ifelse(dim(x)[1] == 1, "split", "splits"),
      if(any(trivial)) paste0('(', sum(trivial), ' trivial)'),
      "dividing", nTip,
      ifelse(is.null(tipLabels), "unlabelled tips.",
             paste("tips,", tipLabels[1], "..", tipLabels[nTip]))
      )
  if (details) {
    splitNames <- rownames(x)
    if (!is.null(splitNames)) {

      nameLengths <- vapply(splitNames, nchar, 0)
      namePads <- vapply(nameLengths, function(thisLength)
        paste0(rep.int(' ', max(nameLengths) - thisLength), collapse=''), character(1))
      splitNames <- paste0(splitNames, namePads)
    } else {
      splitNames <- character(length(x))
      nameLengths = 0L
    }
    cat("\n ", paste0(rep.int(' ', max(nameLengths)), collapse = ''),
        paste0(rep_len(c(1:9, ' '), nTip), collapse = ''))

    for (i in seq_len(dim(x)[1])) {
      split <- x[i, , drop = FALSE]
      cat("\n", splitNames[i], '',
          paste(ifelse(as.logical(rawToBits(split)[seq_len(nTip)]), '*', '.'),
                collapse = ''))
    }
  }
}


#' @family Splits operations
#' @export
summary.Splits <- function(object, ...) {
  print(object, details = TRUE, ...)
  nTip <- attr(object, 'nTip')
  if (is.null(attr(object, 'tip.label'))) {
    cat("\n\nTips not labelled.")
  } else {
    cat("\n\n", paste0("Tip ", seq_len(nTip), ": ", attr(object, 'tip.label'),
                     "\t", c(character(4L), '\n')[seq_len(min(nTip, 5L))]))
  }
}

#' @family Splits operations
#' @export
names.Splits <- rownames


#' @family Splits operations
#' @export
as.character.Splits <- function(x, ...) {
  tipLabels <- attr(x, 'tip.label')
  nTip <- attr(x, 'nTip')

  apply(as.logical(x), 1L, function(inSplit) {
    paste0(paste(tipLabels[inSplit], collapse=' '), ' | ',
           paste(tipLabels[!inSplit], collapse=' '))
  })

}

#' @family Splits operations
#' @export
as.phylo.Splits <- function(x, ...) {
  ret <- structure(list(edge = splits_to_edge(x, NTip(x)),
                        Nnode = NA,
                        tip.label = TipLabels(x)),
                   order = 'preorder',
                   class = 'phylo')
  ret[["Nnode"]] <- dim(ret[["edge"]])[1] + 1 - NTip(ret)
  ret
}

#' @family Splits operations
#' @export
names.Splits <- function(x) rownames(x)

#' @family Splits operations
#' @export
c.Splits <- function(...) {
  splits <- list(...)
  nTip <- unique(vapply(splits, attr, 1, 'nTip'))
  if (length(nTip) > 1L) {
    stop("Splits must relate to identical tips.")
  }
  tips <- lapply(splits, attr, 'tip.label')
  if (length(unique(lapply(tips, sort))) > 1L) {
    stop("All splits must bear identical tips")
  }
  tipLabels <- tips[[1]]
  splits <- c(splits[1], lapply(splits[seq_along(splits)[-1]], as.Splits,
                                  tipLabels = tipLabels))

  structure(do.call(rbind, splits),
            nTip = nTip,
            tip.label = tipLabels,
            class = 'Splits')
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

#' @family Splits operations
#' @export
`!.Splits` <- function(x) {
  nTip <- attr(x, 'nTip')
  x[] <- !x[]
  remainder <- (8L - nTip) %% 8L
  if (remainder) {
    lastSplit <- dim(x)[2]
    endMask <- packBits(c(!logical(8L - remainder), logical(remainder)))
    x[, lastSplit] <- x[, lastSplit] & endMask
  }
  x
}


#' @family Splits operations
#' @export
length.Splits <- function(x) nrow(x)

#' @family Splits operations
#' @importFrom stats setNames
#' @export
duplicated.Splits <- function(x, incomparables = FALSE, ...) {
  dupX <- !x
  useOrig <- x[, 1] < dupX[, 1]
  dupX[useOrig, ] <- x[useOrig, ]
  
  if (dim(dupX)[2] == 1) {
    setNames(duplicated.default(dupX, incomparables = incomparables, ...),
             names(dupX))
  } else {
    duplicated.array(dupX, MARGIN = 1, incomparables = incomparables, ...)
  }
}

#' @family Splits operations
#' @export
unique.Splits <- function(x, incomparables = FALSE, ...) {
  at <- attributes(x)
  dups <- duplicated(x, incomparables, ...)
  x <- x[!dups, , drop = FALSE]
  at[["dim"]] <- dim(x)
  at[["dimnames"]] <- dimnames(x)
  attributes(x) <- at
  x
}

#' Split matching
#'
#' `match()` returns a vector of the positions of (first) matches of splits in
#' its first argument in its second.
#' `%in%` is a more intuitive interface as a binary operator, which returns
#' a logical vector indicating whether there is a match or not for each
#' split in its left operand.
#'
#' `in.Splits()` is an alias for `%in%`, included for backwards compatibility.
#' It will be deprecated in a future release.
#'
#' @param x,table Object of class `Splits`.
#' @param \dots Specify `nomatch =` to provide an integer value that will be
#' used in place of `NA` in the case where no match is found.
# @param incomparables A vector of values that cannot be matched. Any value in
# `x` matching a value in this vector is assigned the `nomatch` value.
# For historical reasons, `FALSE` is equivalent to `NULL`.
#'
#' @return `match()` returns an integer vector specifying the position in
#' `table` that matches each element in `x`, or `nomatch` if no match is found.
#'
#' @examples
#' splits1 <- as.Splits(BalancedTree(7))
#' splits2 <- as.Splits(PectinateTree(7))
#'
#' match(splits1, splits2)
# Turns base functions into S3 generics that can handle `Splits` objects
# (as well as `integer64`s).  This follows equivalent functions in the
# '\pkg{bit64}' package.
#
#' @seealso Corresponding base functions are documented in
#' [`match()`][base::match].
#'
#' @export
#' @keywords methods
# Following https://github.com/cran/bit64/blob/master/R/patch64.R
"match" <- if (!exists("match.default")) {
  function(x, table, ...) UseMethod("match")
} else {
  match
}

#' @method match default
#' @export
"match.default" <- if (!exists("match.default")) {
  function(x, table, ...) base::"match"(x, table, ...)
} else {
  match.default
}

#' @rdname match
#' @family Splits operations
#' @method match Splits
#' @export
match.Splits <- function(x, table, ...) {
  nomatch <- as.integer(c(...)['nomatch'])
  if (length(nomatch) != 1L) {
    nomatch <- NA_integer_
  }

  vapply(seq_along(x), function(i) {
    ret <- which(table %in% x[[i]])
    if (length(ret) == 0) ret <- nomatch
    ret
  }, integer(1))
}

#' @rdname match
#' @export
match.list <- function(x, table, ...) {
  if (inherits(x, 'Splits')) {
    match.Splits(x, table, ...)
  } else {
    NextMethod()
  }
}

#' @rdname match
#' @export
#' @keywords methods
`%in%` <- if (!exists("%in%.default")) {
  function(x, table) UseMethod("%in%")
} else {
  `%in%`
}

#' @method %in% default
#' @export
"%in%.default" <- if (!exists("%in%.default")) {
  function(x, table) base::"%in%"(x, table)
} else {
  `%in%.default`
}

#' @rdname match
#'
#' @return `%in%` returns a logical vector specifying which of the splits in
#' `x` are present in `table`.
#'
#' @examples
#' splits1 %in% splits2
#'
#' @method %in% Splits
#' @export
`%in%.Splits` <- function(x, table) {
  duplicated(c(x, table), fromLast = TRUE)[seq_along(x)]
}

#' @rdname match
#' @export
in.Splits <- `%in%.Splits`

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
  nTip <- attr(x, 'nTip')
  if (!is.numeric(pole)) {
    pole <- match(pole, attr(x, 'tip.label'))
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

