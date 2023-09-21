#' @examples
#' # An example forest of 100 trees, some identical
#' forest <- as.phylo(c(1, rep(10, 79), rep(100, 15), rep(1000, 5)), nTip = 9)
#'
#' # Generate an 80% consensus tree
#' cons <- ape::consensus(forest, p = 0.8)
#' plot(cons)
#' 
#' # Calculate split frequencies
#' splitFreqs <- SplitFrequency(cons, forest)
#'
#' # Optionally, colour edges by corresponding frequency.
#' # Note that not all edges are associated with a unique split
#' # (and two root edges may be associated with one split - not handled here)
#' edgeSupport <- rep(1, nrow(cons$edge)) # Initialize trivial splits to 1
#' childNode <- cons$edge[, 2]
#' edgeSupport[match(names(splitFreqs), childNode)] <- splitFreqs / 100
#' 
#' plot(cons, edge.col = SupportColour(edgeSupport), edge.width = 3)
#' 
#' # Annotate nodes by frequency 
#' LabelSplits(cons, splitFreqs, unit = "%",
#'             col = SupportColor(splitFreqs / 100),
#'             frame = "none", pos = 3L)
#'
