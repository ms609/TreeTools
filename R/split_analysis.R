#' Tips contained within splits
#'
#' `TipsInSplits()` specifies the number of tips that occur within each
#' bipartition split in a `Splits` object.
#'
#' @param splits Object of class `Splits` or `phylo`.
#' @param keep.names Logical specifying whether to include the names of `splits`
#' in the output.
#' @param smallest Logical; if `TRUE`, return the number of leaves in the
#' smaller bipartition.
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
TipsInSplits <- function(splits, keep.names = TRUE, smallest = FALSE, ...) {
  UseMethod('TipsInSplits')
}

#' @rdname TipsInSplits
#' @export
TipsInSplits.Splits <- function(splits, keep.names = TRUE, smallest = FALSE,
                                ...) {
  ret <- tips_in_splits(splits)
  if (smallest) {
    nTip <- NTip(splits)
    big <- ret > nTip / 2
    ret[big] <- nTip - ret[big]
  }
  if (keep.names) names(ret) <- names(splits)
  ret
}

#' @rdname TipsInSplits
#' @export
TipsInSplits.phylo <- function(splits, keep.names = TRUE, smallest = FALSE,
                               ...) {
  TipsInSplits(as.Splits(splits, ...),
               keep.names = keep.names,
               smallest = smallest)
}

#' @rdname TipsInSplits
#' @return `SplitImbalance()` returns a named vector of integers, specifying the
#' number of leaves within a split that are not 'balanced' by a leaf outside it;
#' i.e. a split that divides leaves evenly has an imbalance of zero; one that
#' splits two tips from ten has an imbalance of 10 - 2 = 8.
#' @export
SplitImbalance <- function(splits, keep.names = TRUE, ...) UseMethod('SplitImbalance')

#' @rdname TipsInSplits
#' @export
SplitImbalance.Splits <- function(splits, keep.names = TRUE, ...) {
  nTip <- NTip(splits)
  inSplit <- TipsInSplits(splits, keep.names = TRUE, ...)
  outSplit <- nTip - inSplit

  # Return:
  abs(inSplit - outSplit)
}

#' @rdname TipsInSplits
#' @export
SplitImbalance.phylo <- SplitImbalance.Splits
