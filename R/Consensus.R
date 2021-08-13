#' Consensus trees
#'
#' Calculates the consensus of a set of trees.
#'
#' @param trees List of trees, optionally of class `multiPhylo`.
#' @param p Proportion of trees that must contain a split for it to be reported
#' in the consensus.  `p = 0.5` gives the majority-rule consensus; `p = 1` (the
#' default) gives the strict consensus.
#' @param check.labels Logical specifying whether to check that all trees have
#' identical labels.  Defaults to `TRUE`, which is slower.
#'
#' @return `Consensus()` returns an object of class `phylo`, rooted as in the
#' first entry of `trees`.
#' @examples
#' Consensus(as.phylo(0:2, 8))
#' @template MRS
#' @family consensus tree functions
#' @export
Consensus <- function (trees, p = 1, check.labels = TRUE) {
  repeat{
    nTip <- NTip(trees)
    if (length(unique(nTip)) > 1) {
      warning("Tree sizes differ; removing leaves not in smallest.")
      trees[] <- lapply(trees, KeepTip, trees[[which.min(nTip)]]$tip.label)
    } else {
      nTip <- nTip[1]
      break
    }
  }
  if (nTip < 4L) {
    return(trees[[1]])
  }
  if (check.labels) {
    trees <- RenumberTips(trees, trees[[1]])
  }
  if (p < 0.5 || p > 1) {
    stop("`p` must be between 0.5 and 1.")
  }
  splits <- as.Splits(consensus_tree(trees, p),
                      tipLabels = TipLabels(trees[[1]]))
  tree1 <- Preorder(trees[[1]])
  edg <- tree1$edge
  root <- edg[DescendantEdges(1, edg[, 1], edg[, 2]), 2]
  root <- root[root <= NTip(tree1)]

  # Return:
  RootTree(as.phylo(splits), root)
}
