#' Paths present in reduced tree
#' 
#' Lists which paths present in a master tree are present when leaves are
#' dropped.
#' 
#' @param paths `data.frame` of paths in master tree, perhaps generated using
#' [`PathLengths()`].
#' @param keptVerts Logical specifying whether each entry is retained in the
#' reduced tree, perhaps generated using [`KeptVerts()`].
#' @param all Logical: if `TRUE`, return all paths that occur in the reduced
#' tree; if `FALSE`, return only those paths that correspond to a single edge.
#' that correspond to edges in the reduced tree.
#' Ignored if `paths` is a matrix.
#' @return `KeptPaths()` returns a logical vector specifying whether each path
#' in `paths` occurs when `keptVerts` vertices are retained.
#' @family tree manipulation
#' @examples 
#' master <- BalancedTree(9)
#' paths <- PathLengths(master)
#' keptTips <- c(1, 5, 7, 9)
#' keptVerts <- KeptVerts(master, keptTips)
#' KeptPaths(paths, keptVerts)
#' paths[KeptPaths(paths, keptVerts, all = FALSE), ]
#' @export
#' @template MRS
KeptPaths <- function(paths, keptVerts, all = TRUE) UseMethod("KeptPaths")

#' @rdname KeptPaths
#' @export
KeptPaths.data.frame <- function(paths, keptVerts, all = TRUE) {
  kept <- apply(matrix(as.integer(unlist(paths[, 1:2])) %in% which(keptVerts),
                       ncol = 2L), 1, all)
  
  if (!all) {
    kept[kept][duplicated(paths[kept, "end"], fromLast = TRUE)] <- FALSE
  }
  # Return:
  kept
}

#' @rdname KeptPaths
#' @export
KeptPaths.matrix <- function(paths, keptVerts, all = TRUE) {
  kept <- array(FALSE, dim = dim(paths))
  kept[keptVerts, keptVerts] <- TRUE
  kept
}
