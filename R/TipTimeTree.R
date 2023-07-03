#' Display time-calibrated tree using tip information only
#' 
#' `TipTimeTree()` plots a phylogenetic tree against time using an
#' _ad hoc_ approach based on dates associated with the leaves.
#'  
#' Nodes are dated to the youngest possible value, plus an additional "buffer"
#' (specified with `minEdge`) to ensure that branching order is readable.
#'  
#' @template treeParam
#' @param tipAges Numeric vector specifying the age (in units-of-time ago)
#' associated with each tip in `tree$tip.label` in turn.
#' Older ages signify earlier tips.
#' @param minEdge Minimum length of edge to allow (in units-of-time)
#' 
#' @return `TipTimedTree()` returns a tree with edge lengths set based on the
#' ages of each tip.
#' 
#' @examples
#' tree <- BalancedTree(6)
#' plot(TipTimedTree(tree, tipAge = 1:6, minEdge = 2))
#' @family utility functions
#' @family tree manipulation
#' @export
TipTimedTree <- function(tree, tipAge, minEdge = 1) {
  po <- Preorder(tree)
  edge <- po$edge
  nEdge <- dim(edge)[1]
  if (nEdge < 2) {
    warning("`tree` must contain multiple edges")
    tree[["edge.length"]] <- rep_len(minEdge, nEdge)
    return(tree)
  }
  nTip <- NTip(tree)
  if (nTip != length(tipAge)) {
    stop("`tipAge` must list one age per leaf in `tree`")
  }
  nNode <- nEdge - nTip
  
  age <- c(tipAge, rep_len(-Inf, max(edge) - nTip))
  
  for (i in nEdge:1) {
    age[edge[i, 1]] <- max(age[edge[i, 2]] + minEdge, age[edge[i, 1]])
  }
  tree$edge.length <- apply(edge, 1, function(i) age[i[1]] - age[i[2]])
  # Return:
  tree
}
