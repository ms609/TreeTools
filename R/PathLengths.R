#' Calculate length of paths between each pair of vertices within tree
#' 
#' Given a weighted rooted tree `tree`, `PathLengths()` returns the distance
#' from each vertex to each of its descendant vertices.
#' 
#' @param tree Original tree of class `phylo`, in [`Preorder`].
#' @param fullMatrix Logical specifying return format; see "value" section`.
#' @param use.na Logical specifying whether to set non-existent paths to `NA`,
#' or to leave uninitialized.  Set to `FALSE` to maximize performance.
#' @return If `fullMatrix = TRUE`, `PathLengths()` returns a square matrix in
#' which entry `[i, j]` denotes the distance from internal node `i` to the
#' descendant vertex `j`. 
#' Vertex pairs without a continuous directed path are denoted `NA` if `use.na`
#' is `TRUE`.
#' If `fullMatrix = FALSE`, `PathLengths()` returns a `data.frame` with three
#' columns: `start` lists the deepest node in each path (i.e. that closest
#' to the root); `end` lists the shallowest node (i.e. that closest to a leaf);
#' `length` lists the total length of that path.
#' @examples
#' tree <- rtree(6)
#' plot(tree)
#' add.scale.bar()
#' nodelabels()
#' tiplabels()
#' PathLengths(tree)
#' @template MRS
#' @export
#' @family tree properties
PathLengths <- function(tree, fullMatrix = FALSE, use.na = TRUE) {
  if (!inherits(tree, "phylo")) {
    stop("`tree` must be an object of class `phylo`.")
  }
  weights <- tree[["edge.length"]]
  if (is.null(weights)) {
    weights <- rep_len(1L, dim(tree[["edge"]])[1])
  }
  mat <- path_lengths(tree[["edge"]], weights,
                      init_nas = isFALSE(fullMatrix) || isTRUE(use.na))
  
  # Return:
  if (fullMatrix) {
    mat
  } else {
    paths <- !is.na(mat)
    data.frame(start = row(paths)[paths],
               end = col(paths)[paths],
               length = mat[paths])
  }
}