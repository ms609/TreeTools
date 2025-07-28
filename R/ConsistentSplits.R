#' Are splits consistent?
#' 
#' @param needle Splits object containing the single split to evaluate
#' @param haystack Splits object, or list thereof, containing the splits to
#' compare against `needle`.
#' @examples
#' splits1 <- as.Splits(BalancedTree(8))
#' splits2 <- as.Splits(PectinateTree(8))
#' summary(splits1[[4]])
#' SplitConsistent(splits1[[4]], splits2)
#' 
#' @template MRS
#' @export
SplitConsistent <- function(needle, haystack) {
  if (inherits(haystack, "Splits")) {
    haystack <- list(haystack)
  }
  if (!is.list(haystack)) {
    stop("`haystack` must be a Splits object, or a list thereof.")
  }
  split_consistent(needle, haystack)
}
