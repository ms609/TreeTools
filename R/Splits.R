#' As splits
#'
#' Converts a phylogenetic tree to an array of bipartition splits.
#'
#' @param tr A tree of class \code{\link[ape:read.tree]{phylo}}.
#' @param tipLabels Character vector specifying sequence in which to order
#' tip labels.  Label order must (currently) match to combine or compare separate
#' `Splits` objects.
#' @param asSplits Logical specifying whether to return a `Splits` object,
#'   or an unannotated two-dimensional array (useful where performance is
#'   paramount).
#' @return Returns an object of class `Splits`, or (if `asSplits = FALSE`) a
#'  two-dimensional array of 32-bit integers, which each bit specifying whether
#'  a tip is a member of the split.  Splits are named according to the node
#'  that defines them.
#'
#' @author Martin R. Smith
#'
#' @examples as.Splits(ape::rtree(6, tip.label=1:6, br=NULL))
#'
#' @importFrom ape reorder.phylo
#' @export
as.Splits <- function (x, tipLabels = NULL, ...) UseMethod('as.Splits')

#' @export
as.Splits.phylo <- function (x, tipLabels = NULL, asSplits = TRUE) {
  if (!is.null(tipLabels)) {
    x <- RenumberTips(x, tipLabels)
  }
  x <- Cladewise(x)
  splits <- cpp_edge_to_splits(x$edge)
  # Return:
  if (asSplits) {
    nEdge <- dim(x$edge)[1]
    nTip <- length(x$tip.label)
    rownames(splits) <- paste0('n', (nTip + 3L):(nEdge + 1L))
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
  if (attr(x, 'tip.label') != tipLabels) {
    stop("Order of tip labels does not match")
  }
}

#' @export
as.logical.Splits <- function (x) {
  nTip = attr(x, 'nTip')
  ret <- t(apply(x, 1, function (split) {
    unlist(.DecodeBinary(split, nTip = nTip, print = FALSE))
  }))
  colnames(ret) <- attr(x, 'tip.label')
  rownames(ret) <- rownames(x)
  ret
}

#' @export
as.Splits.list <- function (x, tipLabels = x[[1]]$tip.label, asSplits = TRUE) {
  if (class(x[[1]]) == 'phylo') {
    lapply(x, as.Splits, tipLabels = tipLabels, asSplits = asSplits)
  } else {
    stop("Unsupported list type.")
  }
}

#' @export
as.Splits.logical <- function (x, tipLabels = colnames(x), ...) {
  powersOf2 <- 2L^(0:31)
  dimX <- dim(x)
  if (is.null(dimX)) {
    nTip <- length(x)
    chunks <- (nTip %/% 32L) + 1L
    remainder <- nTip %% 32L

    if (is.null(tipLabels)) {
      tipLabels <- names(x)
      if (is.null(tipLabels)) {
        tipLabels <- paste0('t', seq_len(nTip))
      }
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
                                  asSplits = TRUE) {
  lapply(x, as.Splits, tipLabels = tipLabels, asSplits = asSplits)
}

#' @export
as.Splits.numeric <- function (x, tipLabels = paste0('t', seq_along(x))) {
  if (any(x >= 2^32)) stop ("Input out of range")
  structure(matrix(x, ncol = 1),
            nTip = length(tipLabels),
            tip.label = tipLabels,
            class = 'Splits')
}

#' @exportClass Splits
#

#' @export
print.Splits <- function (x, details = FALSE, ...) {
  nTip <- attr(x, 'nTip')
  cat(dim(x)[1], "bipartition", ifelse(dim(x)[1] > 1, "splits", "split"),
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

#' @export
summary.Splits <- function (x, ...) {
  print(x, details = TRUE, ...)
  nTip <- attr(x, 'nTip')
  cat("\n\n", paste0("Tip ", seq_len(nTip), ": ", attr(x, 'tip.label'),
                   "\t", c(rep('', 4), '\n')[seq_len(min(nTip, 5L))]))
}

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
#' @export
Ntip.Splits <- function (x, ...) attr(x, 'nTip')

#' @export
IdentifySplits <- function (splits, taxonNames = rownames(splits)) {
  apply(splits, 2, function (x) paste0(
    paste(taxonNames[x], collapse=' '), ' : ',
    paste(taxonNames[!x], collapse=' ')))
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
  c(lapply(n[-nN], .DecodeBinary32, print = print, appendLF = appendLF),
    .DecodeBinary32Last(n[nN], nTip, print = print, appendLF = appendLF))
}

#' @export
c.Splits <- function (...) {
  splits <- list(...)
  nTip <- unique(vapply(splits, attr, 1, 'nTip'))
  if (length(nTip) > 1) {
    stop("Splits must relate to identical tips.")
  }
  tips <- vapply(splits, attr, character(nTip), 'tip.label')
  if (dim(unique(tips, MARGIN = 2))[2] != 1) {
    stop("Order of tip labels must be identical.")
  }

  x <- rbind(...)
  structure(x,
            nTip = nTip,
            tip.label = tips[, 1],
            class='Splits')
}


#' @export
length.Splits <- function (x) nrow(x)

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

#' @export
match.Splits <- function (x, table, nomatch = NA_integer_,
                          incomparables = NULL) {
  ret <- which(in.Splits(x, table, incomparables))
  ret[is.na(ret)] <- nomatch
  ret
}

#' @export
in.Splits <- function (x, table, incomparables = NULL, ...) {
  duplicated(c(x, table), fromLast = TRUE,
             incomparables = incomparables)[seq_along(x), ]
}
