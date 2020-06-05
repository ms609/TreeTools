#' @examples
#' # An example forest of 100 trees, some identical
#' forest <- as.phylo(c(1, rep(10, 79), rep(100, 15), rep(1000, 5)), nTip = 9)
#'
#' # Generate an 80% consensus tree
#' cons <- ape::consensus(forest, p = 0.8)
#' plot(cons)
#'
#' splitFreqs <- SplitFrequency(cons, forest)
#' LabelSplits(cons, splitFreqs, unit = '%',
#'             col = SupportColor(splitFreqs / 100),
#'             frame = 'none', pos = 3L)
#'
