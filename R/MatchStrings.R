#' Check for mismatch between character vectors
#' 
#' Checks that entries in one character vector occur in another, suggesting
#' corrections for mismatched elements.
#' 
#' @param x,table Character vectors, in which all elements of `x` are expected
#' to occur in `table`.
#' @param Fail Function to call if a mismatch is found.
#' @inheritParams agrep
#' @param \dots Additional arguments to \code{\link[base:agrep]{agrep()}}.
MatchStrings <- function(x, table, Fail = stop, max.distance = 0.5, ...) {
  matches <- match(x, table)
  missing <- is.na(matches)
  if (any(missing)) {
    nearMiss <- unlist(lapply(x[missing], agrep, table,
                              max.distance = max.distance, ...),
                       use.names = FALSE, recursive = FALSE)
    message <- paste0("Could not find '", paste(x[missing], collapse = "', '"), 
                      "' in ", deparse(substitute(table)), ".  ",
                      if (length(nearMiss)) {
                        paste0("Did you mean '", 
                               paste(table[nearMiss], collapse = "', '"), "'?")
                      })
    Fail(message)
  }
}