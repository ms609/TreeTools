#' As splits
#'
#' Converts a phylogenetic tree to an array of bipartition splits.
#'
#' @param tr A tree of class \code{\link[ape:read.tree]{phylo}}.
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
#' @examples Tree2Splits(ape::rtree(6, tip.label=1:6, br=NULL))
#'
#' @importFrom ape reorder.phylo
#' @export
as.Splits <- function (x, ...) UseMethod('as.Splits')

as.Splits.phylo <- function (x, asSplits = TRUE) {
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
  cat(dim(x)[1], "bipartition splits dividing", nTip, "tips.")
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
      lapply(split[-length(split)], .DecodeBinary32)
      .DecodeBinary32(split[length(split)],
                      stopAt = ifelse(nTip %% 32, nTip %% 32, 32))
    }
  }
}

#' @export
summary.Splits <- function (x, ...) {
  print(x, details = TRUE, ...)
}


#' @export
names.Splits <- function (x) rownames(x)

#' @keywords internal
#' @export
.DecodeBinary32 <- function (n, stopAt = 32L, appendLF = FALSE) {
  for (i in seq_len(stopAt)) {
    cat(ifelse(n %% 2, '*', '.'))
    n <- n %/% 2
  }
  if (appendLF) cat("\n")
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
unique.Splits <- function (x, incomparables = FALSE) {
  at <- attributes(x)
  nTip <- at$nTip
  notX <- (2^32 - 1) - x
  notX[, 1] <- (2 ^ nTip - 1) - x[, 1]
  dupX <- x
  useNot <- notX[, 1] < x[, 1]
  dupX[useNot, ] <- notX[useNot, ]
  dups <- duplicated.array(dupX, MARGIN = 1, incomparables = incomparables)
  x <- x[!dups, , drop = FALSE]
  at$dim <- dim(x)
  at$dimnames <- dimnames(x)
  attributes(x) <- at
  x
}


#' Drop Single Splits
#'
#' Removes splits that pertain only to a single taxon from a splits object.
#'
#' Bipartition splits are divisions, implied by each edge or node of an unrooted
#' tree topology, that divide the taxa into two groups (one of which is a clade).
#'
#' By default, a list of splits will include those that separate a single taxon
#' (a leaf) from all others.  Such splits are, by definition, present in all
#' trees that contain that taxon; they are not of interest when comparing trees.
#' This function removes such splits from a list of bipartitions.
#'
#' @param split A matrix in which each column corresponds to a bipartition split
#'
#' @return The input matrix, with any columns that separate only a single pendant
#'  tip removed.
#'
#' @author Martin R. Smith
#'
#' @export
DropSingleSplits <- function (split) {
  split[, colSums(split) > 1 & colSums(!split) > 1, drop=FALSE]
}
