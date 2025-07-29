#' Identify consistent / conflicting splits
#' 
#' Check whether a series of splits are consistent with or contradict a focal
#' split.
#' 
#' @param needle Splits object containing the single split to evaluate
#' @param haystack Splits object, or list thereof, containing the splits to
#' compare against `needle`.
#' @returns `SplitConsistent()` returns a list of logical vectors.
#' Each list item corresponds to an entry in `haystack`, reporting whether each
#' split is consistent with (`TRUE`) or in conflict with (`FALSE`) `needle`.
#' `SplitConflicts()` returns the inverse.
#' @examples
#' splits1 <- as.Splits(BalancedTree(8))
#' splits2 <- as.Splits(PectinateTree(8))
#' summary(splits1[[4]])
#' SplitConsistent(splits1[[4]], splits2)
#' @template MRS
#' @family split manipulation functions
#' @export
SplitConsistent <- function(needle, haystack) {
  if (inherits(haystack, "Splits")) {
    haystack <- list(haystack)
  }
  if (!is.list(haystack)) {
    stop("`haystack` must be a Splits object, or a list thereof.")
  }
  split_consistent(needle, haystack, FALSE)
}

#' @rdname SplitConsistent
#' @export
SplitConflicts <- function(needle, haystack) {
  if (inherits(haystack, "Splits")) {
    haystack <- list(haystack)
  }
  if (!is.list(haystack)) {
    stop("`haystack` must be a Splits object, or a list thereof.")
  }
  split_consistent(needle, haystack, TRUE)
}
