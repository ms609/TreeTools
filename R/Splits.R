#' As splits
#'
#' Converts a phylogenetic tree to an array of bipartition splits.
#'
#' @param x Object to convert into splits: perhaps a tree of class
#'  \code{\link[ape:read.tree]{phylo}}.
#' @param tipLabels Character vector specifying sequence in which to order
#' tip labels.  Label order must (currently) match to combine or compare separate
#' `Splits` objects.
#' @param \dots Presently unused.
#' @return Returns an object of class `Splits`, or (if `asSplits = FALSE`) a
#'  two-dimensional array of 32-bit integers, which each bit specifying whether
#'  a tip is a member of the split.  Splits are named according to the node
#'  that defines them.
#'
#' @author Martin R. Smith
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
#' @export
as.Splits <- function (x, tipLabels = NULL, ...) UseMethod('as.Splits')

#' @describeIn as.Splits Convert object of class `phylo` to `Splits`.
#' @param asSplits Logical specifying whether to return a `Splits` object,
#'   or an unannotated two-dimensional array (useful where performance is
#'   paramount).
#' @export
as.Splits.phylo <- function (x, tipLabels = NULL, asSplits = TRUE, ...) {
  if (!is.null(tipLabels)) {
    x <- RenumberTips(x, .TipLabels(tipLabels))
  }
  x <- Cladewise(x)
  splits <- cpp_edge_to_splits(x$edge)
  nSplits <- dim(splits)[1]
  # Return:
  if (asSplits) {
    nEdge <- dim(x$edge)[1]
    nTip <- length(x$tip.label)
    if (nSplits > 0) rownames(splits) <- paste0('n', nTip + 2L + seq_len(nSplits))
    structure(splits,
              nTip = nTip,
              tip.label = x$tip.label,
              class = 'Splits')
  }
  else {
    splits
  }
}

#' @export
as.Splits.Splits <- function (x, tipLabels = NULL, ...) {
  tipLabels <- .TipLabels(tipLabels)
  oldLabels <- attr(x, 'tip.label')
  nTip <- attr(x, 'nTip')
  if (is.null(oldLabels)) {
    if (length(tipLabels) == nTip) {
      attr(x, 'tip.label') <- tipLabels
      x
    } else {
      stop (length(tipLabels), " labels provided; expecting ", nTip)
    }
  }
  if (!identical(oldLabels, tipLabels)) {
    if (all(oldLabels %in% tipLabels) && all(tipLabels %in% oldLabels)) {
      ret <- as.Splits(t(apply(x, 1, .DecodeBinary, nTip = nTip)
                         [match(tipLabels, oldLabels), ]))
      attr(ret, 'tip.label') <- tipLabels
      ret
    } else {
      stop ("Old and new labels must match")
    }
  } else {
    x
  }
}

#' @export
as.logical.Splits <- function (x, tipLabels = NULL, ...) {
  nTip = attr(x, 'nTip')
  ret <- t(apply(x, 1, function (split) {
    unlist(.DecodeBinary(split, nTip = nTip, print = FALSE))
  }))
  colnames(ret) <- attr(x, 'tip.label')
  rownames(ret) <- rownames(x)
  ret
}

#' @export
as.Splits.list <- function (x, tipLabels = x[[1]]$tip.label, asSplits = TRUE, ...) {
  if (class(x[[1]]) == 'phylo') {
    lapply(x, as.Splits, tipLabels = tipLabels, asSplits = asSplits)
  } else {
    stop("Unsupported list type.")
  }
}

#' @export
as.Splits.logical <- function (x, tipLabels = NULL, ...) {
  powersOf2 <- 2L^(0:31)
  dimX <- dim(x)
  if (is.null(dimX)) {
    nTip <- length(x)
    chunks <- (nTip %/% 32L) + 1L
    remainder <- nTip %% 32L

    if (is.null(tipLabels)) {
      tipLabels <- .TipLabels(x)
      if (is.null(tipLabels)) {
        tipLabels <- paste0('t', seq_len(nTip))
      }
    } else {
      tipLabels <- .TipLabels(tipLabels)
    }

    structure(matrix(vapply(seq_len(chunks) - 1L, function (i) {
      chunk <- seq_len(
        if (i + 1L == chunks && remainder != 0L) remainder else 32L
      )
      sum(powersOf2[chunk][x[i * 32L + chunk]])
    }, double(1)), nrow = 1L),
    nTip = nTip,
    tip.label = tipLabels,
    class = 'Splits')
  } else {
    nTip <- dimX[2]
    chunks <- (nTip %/% 32L) + 1L
    remainder <- nTip %% 32L

    structure(vapply(seq_len(chunks) - 1L, function (i) {
      chunkSeq <- seq_len(
        if (i + 1L == chunks && remainder != 0L) remainder else 32L
      )
      chunk <- i * 32L + chunkSeq
      apply(x[, chunk, drop = FALSE], 1L,
            function (twos) sum(powersOf2[chunkSeq][twos]))
    }, double(dimX[1])),
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

#' @export
as.Splits.numeric <- function (x, tipLabels = paste0('t', seq_along(x)), ...) {
  if (any(x >= 2^32)) stop ("Input out of range")
  structure(matrix(x, ncol = 1),
            nTip = length(tipLabels),
            tip.label = tipLabels,
            class = 'Splits')
}

#' @exportClass Splits
#

#' @family Splits operations
#' @export
print.Splits <- function (x, details = FALSE, ...) {
  nTip <- attr(x, 'nTip')
  cat(dim(x)[1], "bipartition", ifelse(dim(x)[1] == 1, "split", "splits"),
      "dividing", nTip, "tips.")
  if (details) {
    #cat(paste0("\n  Tip ", seq_len(nTip), ': ', attr(x, 'tip.label')))
    #for (i in rev(seq_len(log10(nTip) + 1L)) - 1L) {
    #  cat("\n", vapply(seq_len(1L + (nTip %/% 10^i)) - 1, function (x) {
    #    rep(ifelse(x == 0, ' ', as.character(x %% 10^i)), 10^i)}
    #    , character(10^i)))
    #}

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
  cat("\n\n", paste0("Tip ", seq_len(nTip), ": ", attr(object, 'tip.label'),
                   "\t", c(rep('', 4), '\n')[seq_len(min(nTip, 5L))]))
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

#' @importFrom ape Ntip
#' @family Splits operations
#' @export
Ntip.Splits <- function (phy, ...) attr(phy, 'nTip')

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

#' @keywords internal
#' @export
.DecodeBinary32 <- function (n, stopAt = 32L, print = FALSE, appendLF = FALSE) {
  ret <- logical(stopAt)
  for (i in seq_len(stopAt)) {
    ret[i] = as.logical(n %% 2)
    if (print) cat(ifelse(ret[i], '*', '.'))
    n <- n %/% 2
  }
  if (print && appendLF) cat("\n")
  ret
}

#' @keywords internal
#' @export
.DecodeBinary32Last <- function (n, nTip, ...) {
  remainder32 <- nTip %% 32L
  .DecodeBinary32(n, stopAt = ifelse(remainder32, remainder32, 32L), ...)
}

#' @keywords internal
#' @export
.DecodeBinary <- function (n, nTip, print = FALSE, appendLF = FALSE, ...) {
  nN <- length(n)
  c(vapply(n[-nN], .DecodeBinary32, print = print, appendLF = appendLF,
           FUN.VALUE = logical(32L)),
    .DecodeBinary32Last(n[nN], nTip, print = print, appendLF = appendLF))
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

#' @keywords internal
#' @export
.LastBin <- function (n) if (n %% 32L) n %% 32L else 32L

#' @keywords internal
#' @export
.BinSizes <- function (n) c(rep(32L, (n - 1L) %/% 32L), if (n %% 32L) n %% 32L else 32L)


#' @family Splits operations
#' @export
`!.Splits` <- function (x) {
  dimX <- dim(x)
  matrix(2L ^ .BinSizes(Ntip(x)) - 1L, dimX[1], dimX[2], byrow=TRUE) - x
}



#' @family Splits operations
#' @export
length.Splits <- function (x) nrow(x)

#' @family Splits operations
#' @export
duplicated.Splits <- function (x, incomparables = FALSE, ...) {
  nTip <- attr(x, 'nTip')
  dimX <- dim(x)

  dupX <- matrix(nrow = dimX[1], ncol = dimX[2])
  dupX[, 1] <- (2 ^ nTip - 1) - x[, 1]
  dupX[, -1] <- (2^32 - 1) - x[, -1]

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
#' @return A logical vector specifing which of the splits in `x` are present
#' in `table`.
#'
#' @author Martin R. Smith
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
             incomparables = incomparables)[seq_along(x), ]
}


#' @keywords internal
#' @export
.TipLabels <- function (x) UseMethod('.TipLabels')

#' @keywords internal
#' @export
.TipLabels.default <- function (x) {
  if (is.null(names(x))) {
    if (any(duplicated(x))) {
      NULL
    } else {
      x
    }
  } else {
    names(x)
  }
}

#' @keywords internal
#' @export
.TipLabels.phylo <- function (x) x$tip.label

#' @keywords internal
#' @export
.TipLabels.list <- function (x) {
  .TipLabels(x[[1]])
}

#' @keywords internal
#' @export
.TipLabels.matrix <- function (x) colnames(x)

#' @keywords internal
#' @export
.TipLabels.multiPhylo <- function (x) {
  .TipLabels(x[[1]])
}

#' @keywords internal
#' @export
.TipLabels.Splits <- function (x) attr(x, 'tip.label')


#' @keywords internal
#' @export
.TipLabels.numeric <- function (x) NextMethod('.TipLabels', as.character(x))
