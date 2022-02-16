#' Identify vertices retained when leaves are dropped
#' 
#' @param tree Original tree of class `phylo`, in [`Preorder`].
#' @param keptTips Logical vector stating whether each leaf should be retained.
#' Sequence corresponds to `tree[["tip.label"]]`.
#' @examples 
#' master <- BalancedTree(12)
#' master <- Preorder(master) # Nodes must be listed in Preorder sequence
#' plot(master)
#' nodelabels()
#' 
#' allTips <- master[["tip.label"]]
#' keptTips <- sample(allTips, 8)
#' plot(KeepTip(master, keptTips))
#' kept <- KeptVerts(master, allTips %in% keptTips)
#' 
#' map <- which(kept)
#' # Node `i` in the reduced tree corresponds to node `map[i]` in the original.
#' @template MRS
#' @export
#' @family tree manipulation
KeptVerts <- function(tree, keptTips) UseMethod("KeptVerts")

#' @rdname KeptVerts
#' @export
KeptVerts.phylo <- function(tree, keptTips) {
  order <- attr(tree, "order")
  if (is.null(order) || order != "preorder") {
    stop("`tree` must be in preorder; try `Preorder(tree)`")
  }
  KeptVerts(tree[["edge"]], keptTips)
}

#' @rdname KeptVerts
#' @export
KeptVerts.numeric <- function(tree, keptTips) {
  if (!is.logical(keptTips)) {
    stop("`keptTips` must be a logical vector")
  }
  dims <- dim(tree)
  if (is.null(dims) || dims[2] != 2L) {
    stop("`tree` must be the numeric edge matrix of a `phylo` object")
  }
  ret <- kept_vertices(tree, keptTips)
  # Return:
  ret[-1] > 1L
}
