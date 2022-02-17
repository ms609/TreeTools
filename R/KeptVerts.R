#' Identify vertices retained when leaves are dropped
#' 
#' @param tree Original tree of class `phylo`, in [`Preorder`].
#' @param keptTips Either: 
#' - a logical vector stating whether each leaf should be retained, in a
#' sequence corresponding to `tree[["tip.label"]]`; or
#' - a character vector listing the leaf labels to retain; or
#' - a numeric vector listing the indices of leaves to retain.
#' @param tipLabels Optional character vector naming the leaves of `tree`,
#' if `keptTips` is not logical.  Inferred from `tree` if unspecified.
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
KeptVerts <- function(tree, keptTips, tipLabels = TipLabels(tree)) {
  UseMethod("KeptVerts")
}

#' @rdname KeptVerts
#' @export
KeptVerts.phylo <- function(tree, keptTips, tipLabels = TipLabels(tree)) {
  order <- attr(tree, "order")
  if (is.null(order) || order != "preorder") {
    stop("`tree` must be in preorder; try `Preorder(tree)`")
  }
  KeptVerts(tree[["edge"]], keptTips, tipLabels = tipLabels)
}

#' @rdname KeptVerts
#' @export
KeptVerts.numeric <- function(tree, keptTips, tipLabels = TipLabels(tree)) {
  switch(mode(keptTips),
         "numeric" = {
           x <- logical(length(tipLabels))
           x[keptTips] <- TRUE
           keptTips <- x
         },
         "character" = {
           keptTips <- tipLabels %in% keptTips
         },
         "logical" = {},
         stop("Unrecognized format for `keptTips`")
  )
  dims <- dim(tree)
  if (is.null(dims) || dims[2] != 2L) {
    stop("`tree` must be the numeric edge matrix of a `phylo` object")
  }
  ret <- kept_vertices(tree, keptTips)
  # Return:
  ret[-1] > 1L
}
