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
#' KeptVerts(master, allTips %in% keptTips)
#' @template MRS
#' @export
#' @family tree manipulation
KeptVerts <- function(tree, keptTips) UseMethod("KeptVerts")

#' @rdname KeptVerts
#' @export
KeptVerts.phylo <- function(tree, keptTips) {
  KeptVerts(tree[["edge"]], keptTips)
}

#' @rdname KeptVerts
#' @export
KeptVerts.matrix <- function(tree, keptTips) {
  ret <- kept_vertices(tree, keptTips)
  # Return:
  ret[-1] > 1L
}
