#' Paths present in reduced tree
#' 
#' Lists which paths present in a master tree are present when leaves are
#' dropped.
#' 
#' @param paths `data.frame` of paths in master tree, perhaps generated using
#' [`PathLengths()`].
#' @param keptVerts Logical specifying whether each entry is retained in the
#' reduced tree, perhaps generated using [`KeptVerts()`].
#' @return `KeptPaths()` returns a logical vector specifying whether each path
#' in `paths` occurs when `keptVerts` vertices are retained.
#' @family tree manipulation
#' @export
#' @template MRS
KeptPaths <- function(paths, keptVerts) UseMethod("KeptPaths")

#' @rdname KeptPaths
#' @export
KeptPaths.data.frame <- function(paths, keptVerts) {
  apply(matrix(as.integer(unlist(paths[, 1:2])) %in% which(keptVerts),
               ncol = 2L), 1, all)
}

#' @rdname KeptPaths
#' @export
KeptPaths.matrix <- function(paths, keptVerts) {
  kept <- array(FALSE, dim = dim(paths))
  kept[keptVerts, keptVerts] <- TRUE
  kept
}
