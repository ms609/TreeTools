#' Tips contained within splits
#'
#' `TipsInSplits()` specifies the number of tips that occur within each
#' bipartition split in a `Splits` object.
#'
#' @param splits Object of class `Splits` or `phylo`.
#' @param \dots Additional parameters to pass to `as.Splits()`.
#'
#' @return `TipsInSplits()` returns a named vector of integers, specifying the
#' number of tips contained within each split in `splits`.
#'
#' @examples
#' tree <- PectinateTree(8)
#' splits <- as.Splits(tree)
#' TipsInSplits(splits)
#'
#' plot(tree)
#' LabelSplits(tree, as.character(splits), frame = 'none', pos = 3L, cex = 0.7)
#' LabelSplits(tree, TipsInSplits(splits), unit = ' tips', frame = 'none',
#'             pos = 1L)
#'
#' @family Splits operations
#' @export
TipsInSplits <- function (splits, ...) UseMethod('TipsInSplits')

#' @rdname TipsInSplits
#' @export
TipsInSplits.Splits <- function (splits, ...) {
  dims <- dim(splits)
  ret <- .colSums(as.logical(rawToBits(t(splits))), 8L * dims[2], dims[1])
  names(ret) <- names(splits)
  ret
}

#' @rdname TipsInSplits
#' @export
TipsInSplits.phylo <- function (splits, ...) {
  TipsInSplits(as.Splits(splits, ...))
}

#' @rdname TipsInSplits
#' @return `SplitImbalance()` returns a named vector of integers, specifying the
#' number of leaves within a split that are not 'balanced' by a leaf outside it;
#' i.e. a split that divides leaves evenly has an imbalance of zero; one that
#' splits two tips from ten has an imbalance of 10 - 2 = 8.
#' @export
SplitImbalance <- function (splits, ...) UseMethod('SplitImbalance')

#' @rdname TipsInSplits
#' @export
SplitImbalance.Splits <- function (splits, ...) {
  nTip <- NTip(splits)
  inSplit <- TipsInSplits(splits, ...)
  outSplit <- nTip - inSplit

  # Return:
  abs(inSplit - outSplit)
}

#' @rdname TipsInSplits
#' @export
SplitImbalance.phylo <- SplitImbalance.Splits
