#' As splits
#'
#' Converts a phylogenetic tree to an array of bipartition splits.
#'
#'
#' @param x Object to convert into splits: perhaps a tree of class
#'  \code{\link[ape:read.tree]{phylo}}.
#'  If a logical matrix is provided, each row will be considered as a
#'  separate split.
#' @param tipLabels Character vector specifying sequence in which to order
#' tip labels.  Label order must (currently) match to combine or compare separate
#' `Splits` objects.
#' @param \dots Presently unused.
#' @return Returns an object of class `Splits`, or (if `asSplits = FALSE`) a
#'  two-dimensional array of 32-bit integers, which each bit specifying whether
#'  a tip is a member of the split.  Splits are named according to the node
#'  that defines them.
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
#' match.Splits(splits, moreSplits)
#' in.Splits(moreSplits, splits)
#'
#' @family Splits operations
#' @importFrom ape reorder.phylo
#' @exportClass Splits
#' @export
as.Splits <- function (x, tipLabels = NULL, ...) UseMethod('as.Splits')

#' @rdname as.Splits
#' @param asSplits Logical specifying whether to return a `Splits` object,
#'   or an unannotated two-dimensional array (useful where performance is
#'   paramount).
#' @export
as.Splits.phylo <- function (x, tipLabels = NULL, asSplits = TRUE, ...) {
  if (!is.null(tipLabels)) {
    x <- RenumberTips(x, TipLabels(tipLabels))
  }
  edge <- Postorder(x)$edge
  nTip <- length(x$tip.label)
  splits <- cpp_edge_to_splits(edge, nTip)

  nSplits <- dim(splits)[1]
  # Return:
  if (asSplits) {
    nEdge <- dim(x$edge)[1]
    nTip <- length(x$tip.label)
    structure(splits,
              nTip = nTip,
              tip.label = x$tip.label,
              class = 'Splits')
  }
  else {
    splits
  }
}

#' @rdname as.Splits
#' @export
as.Splits.Splits <- function (x, tipLabels = NULL, ...) {
  if (is.null(tipLabels)) {
    # Nothing needs doing
    # Return:
    x
  } else {
    tipLabels <- TipLabels(tipLabels)
    oldLabels <- attr(x, 'tip.label')
    if (is.null(oldLabels)) {
      nTip <- attr(x, 'nTip')
      if (length(tipLabels) == nTip) {
        attr(x, 'tip.label') <- tipLabels
        x
      } else {
        stop (length(tipLabels), " labels provided; expecting ", nTip)
      }
    }
    if (!identical(oldLabels, tipLabels)) {
      nTip <- attr(x, 'nTip')
      if (length(x) == 0) {
        attr(x, 'tip.label') <- tipLabels
        x
      } else {
        if (all(oldLabels %in% tipLabels) && all(tipLabels %in% oldLabels)) {
          as.Splits(t(apply(x, 1, .DecodeBinary, nTip = nTip)
                             [match(tipLabels, oldLabels), ]),
                           tipLabels = tipLabels)
        } else {
          stop ("Old and new labels must match")
        }
      }
    } else {
      x
    }
  }
}

#' @rdname as.Splits
#' @export
as.Splits.list <- function (x, tipLabels = NULL, asSplits = TRUE, ...) {
  if (inherits(x[[1]], 'phylo')) {
    if (is.null(tipLabels)) {
      tipLabels <- x[[1]]$tip.label
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

#' @rdname as.Splits
#' @export
as.Splits.logical <- function (x, tipLabels = NULL, ...) {
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

    structure(matrix(packBits(c(x, rep(F, (8L - nTip) %% 8))), nrow = 1L),
    nTip = nTip,
    tip.label = tipLabels,
    class = 'Splits')
  } else {
    nTip <- dimX[2]
    if (is.null(tipLabels)) {
      tipLabels <- TipLabels(x)
    }
    if (is.null(tipLabels)) {
      tipLabels <- paste0('t', seq_len(nTip))
    }

    structure(matrix(packBits(t(cbind(x, matrix(F, dimX[1], (8L - nTip) %% 8)))),
                     nrow = dimX[1], byrow=TRUE, dimnames = list(rownames(x), NULL)),
      nTip = nTip,
      tip.label = tipLabels,
      class = 'Splits')
  }
}

#' @export
as.Splits.multiPhylo <- function (x, tipLabels = x[[1]]$tip.label,
                                  asSplits = TRUE, ...) {
  lapply(x, as.Splits, tipLabels = tipLabels, asSplits = asSplits)
}


#' @rdname as.Splits
#' @export
as.logical.Splits <- function (x, tipLabels = NULL, ...) {
  nTip <- attr(x, 'nTip')
  ret <- t(apply(x, 1, function (split) {
    unlist(.DecodeBinary(split, nTip = nTip, print = FALSE))
  }))
  colnames(ret) <- attr(x, 'tip.label')
  rownames(ret) <- rownames(x)
  ret
}

#' @family Splits operations
#' @export
print.Splits <- function (x, details = FALSE, ...) {
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
    names <- rownames(x)
    if (!is.null(names)) {

      nameLengths <- vapply(names, nchar, 0)
      namePads <- vapply(nameLengths, function (thisLength)
        paste0(rep(' ', max(nameLengths) - thisLength), collapse=''), character(1))
      names <- paste0(names, namePads)
    } else {
      names <- rep('', length(x))
      nameLengths = 0L
    }
    cat("\n ", paste0(rep(' ', max(nameLengths)), collapse = ''),
        paste0(rep(c(1:9, ' '), length.out = nTip), collapse=''))

    for (i in seq_len(dim(x)[1])) {
      split <- x[i, , drop = FALSE]
      cat("\n", names[i], ' ')
      .DecodeBinary(split, nTip, print = TRUE)
    }
  }
}


#' @family Splits operations
#' @export
summary.Splits <- function (object, ...) {
  print(object, details = TRUE, ...)
  nTip <- attr(object, 'nTip')
  if (is.null(attr(object, 'tip.label'))) {
    cat("\n\nTips not labelled.")
  } else {
    cat("\n\n", paste0("Tip ", seq_len(nTip), ": ", attr(object, 'tip.label'),
                     "\t", c(rep('', 4), '\n')[seq_len(min(nTip, 5L))]))
  }
}


#' @family Splits operations
#' @export
as.character.Splits <- function (x, ...) {
  tipLabels <- attr(x, 'tip.label')
  nTip <- attr(x, 'nTip')

  apply(x, 1, function (split) {
    inSplit <- unlist(.DecodeBinary(split, nTip))
    paste0(paste(tipLabels[inSplit], collapse=' '), ' | ',
           paste(tipLabels[!inSplit], collapse=' '))

  })
}


#' Number of tips in a phylogenetic tree
#'
#' Extends ape's function [`Ntip`][ape::summary.phylo] to handle objects of
#' class `Splits` and `list`.
#'
#' @param phy Object to count.
#'
#' @return Single integer specifying the number of tips in each object in `phy`.
#'
#' @export
NTip <- function (phy) UseMethod('NTip')

NTip.default <- function (phy) UseMethod('Ntip')

#' @rdname NTip
#' @family Splits operations
#' @export
NTip.Splits <- function (phy) attr(phy, 'nTip')

#' @rdname NTip
#' @export
NTip.list <- function (phy) vapply(phy, NTip, integer(1))

#' @rdname NTip
#' @export
NTip.phylo <- function (phy) length(phy$tip.label)

#' @rdname NTip
#' @export
NTip.multiPhylo <- function (phy) {
  ret <- attr(phy, 'TipLabel')
  if (is.null(ret)) {
    vapply(phy, NTip.phylo, integer(1))
  } else {
    rep(length(ret), length(phy))
  }
}

#' Tips contained within splits
#'
#' `TipsInSplits` specifies the number of tips that occur within each
#' bipartition split in a `Splits` object.
#'
#' @param splits Object of class `Splits`.
#' @param nTip Number of tips in `Splits` object (inferred if not specified).
#'
#' @return A named vector of integers, specifying the number of tips contained
#' within each split in `splits`.
#'
#' @examples
#' splits <- as.Splits(PectinateTree(8))
#' TipsInSplits(splits)
#'
#' @family Splits operations
#' @export
TipsInSplits <- function (splits, nTip = attr(splits, 'nTip')) {
  apply(splits, 1, function (split) sum(unlist(.DecodeBinary(split, nTip = nTip))))
}

#' @export
names.Splits <- function (x) rownames(x)

#' @rdname Decoders
#' @keywords internal
#' @export
.DecodeRaw <- function (n, stopAt = 8L, print = FALSE, appendLF = FALSE) {
  bitMasks <- as.raw(c(1, 2, 4, 8, 16, 32, 64, 128)[seq_len(stopAt)])
  ret <- as.logical(n & bitMasks)
  if (print) cat(paste0(ifelse(ret, '*', '.'), collapse = ''))
  if (print && appendLF) cat("\n")
  ret
}

#' @rdname Decoders
#' @keywords internal
#' @export
.DecodeLastRaw <- function (n, nTip, ...) {
  remainder <- nTip %% 8L
  .DecodeRaw(n, stopAt = ifelse(remainder, remainder, 8L), ...)
}

#' Decode Splits objects
#'
#' Internal functions to decode raw representation of splits.
#'
#' @name Decoders
#' @param n Number to decode.
#' @param stopAt Integer specifying number of tips in partial raw element.
#' @param nTip Integer specifying number of tips in original tree.
#' @param print Logical specifying whether output is to be printed.
#' @param appendLF Logical specifying whether to append line feed character to
#'  output.
#'
#' @return Return a human-readable representation of split content.
#'
#' @keywords internal
#' @export
.DecodeBinary <- function (n, nTip, print = FALSE, appendLF = FALSE, ...) {
  nN <- length(n)
  c(vapply(n[-nN], .DecodeRaw, print = print, appendLF = appendLF,
           FUN.VALUE = logical(8L)),
    .DecodeLastRaw(n[nN], nTip, print = print, appendLF = appendLF))
}

#' @family Splits operations
#' @export
c.Splits <- function (...) {
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
            class='Splits')
}

#' @family Splits operations
#' @export
`+.Splits` <- function (...) c(...)


#' @family Splits operations
#' @export
`-.Splits` <- function (x, y) {
  x[[!in.Splits(x, y)]]
}

#' @family Splits operations
#' @export
`[[.Splits` <- function (x, ..., drop = TRUE) {
  elements <- list(...)
  if (length(elements) == 0L) {
    x[]
  } else if (length(elements) == 1L) {
    ret <- x[elements[[1]], , drop = FALSE]
    at <- attributes(x)
    at$dim<- dim(ret)
    at$dimnames[[1]] <- rownames(x)[elements[[1]]]
    attributes(ret) <- at
    ret
  } else {
    stop ("Too many dimensions specified")
  }
}

#' @family Splits operations
#' @export
`!.Splits` <- function (x) {
  nTip <- attr(x, 'nTip')
  x <- !unclass(x)
  remainder <- (8L - nTip) %% 8L
  if (remainder) {
    lastSplit <- dim(x)[2]
    endMask <- packBits(c(rep(TRUE, 8L - remainder), rep(FALSE, remainder)))
    x[, lastSplit] <- x[, lastSplit] & endMask
  }
  class(x) <- 'Splits'
  x
}


#' @family Splits operations
#' @export
length.Splits <- function (x) nrow(x)

#' @family Splits operations
#' @export
duplicated.Splits <- function (x, incomparables = FALSE, ...) {
  dupX <- !x
  useOrig <- x[, 1] < dupX[, 1]
  dupX[useOrig, ] <- x[useOrig, ]

  duplicated.array(dupX, MARGIN = 1, ...)
}

#' @family Splits operations
#' @export
unique.Splits <- function (x, incomparables = FALSE, ...) {
  at <- attributes(x)
  dups <- duplicated(x, incomparables, ...)
  x <- x[!dups, , drop = FALSE]
  at$dim <- dim(x)
  at$dimnames <- dimnames(x)
  attributes(x) <- at
  x
}

#' Match splits
#'
#' Equivalent of `match` for `Splits` objects.
#'
#' @param x,table Object of class `Splits`.
#' @param nomatch The value to be returned in the case where no match is found.
#' @param incomparables A vector of values that cannot be matched. Any value in
#' `x` matching a value in this vector is assigned the `nomatch` value.
#' For historical reasons, `FALSE` is equivalent to `NULL`.
#'
#' @return An integer vector specifying the position in `table` that matches
#' each element in `x`, or `nomatch` if no match is found.
#'
#' @examples
#' splits1 <- as.Splits(BalancedTree(7))
#' splits2 <- as.Splits(PectinateTree(7))
#'
#' match.Splits(splits1, splits2)
#'
#' @family Splits operations
#' @export
match.Splits <- function (x, table, nomatch = NA_integer_,
                          incomparables = NULL) {
  vapply(seq_along(x), function (i) {
    ret <- which(in.Splits(table, x[[i]], incomparables))
    if (length(ret) == 0) ret <- nomatch
    ret
  }, integer(1))
}

#' Splits in Splits object
#'
#' `in.Splits` is an equivalent to `%in%` that can be applied to objects
#' of class `Splits`.
#'
#' @param x,table Object of class `Splits`.
#' @param incomparables A vector of values that cannot be matched. Any value in
#' `x` matching a value in this vector is assigned the `nomatch` value.
#' For historical reasons, `FALSE` is equivalent to `NULL`.
#'
#' @return A logical vector specifying which of the splits in `x` are present
#' in `table`.
#'
#' @template MRS
#'
#' @examples
#' splits1 <- as.Splits(BalancedTree(7))
#' splits2 <- as.Splits(PectinateTree(7))
#'
#' in.Splits(splits1, splits2)
#'
#' @family Splits operations
#' @export
in.Splits <- function (x, table, incomparables = NULL) {
  duplicated(c(x, table), fromLast = TRUE,
             incomparables = incomparables)[seq_along(x)]
}
