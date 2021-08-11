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
#' @template MRS
#' @family consensus tree functions
#' @export
Consensus <- function (trees, p = 1, check.labels = TRUE) {
  if (check.labels) {
    trees <- RenumberTips(trees, trees[[1]])
  }
  if (p < 0.5 || p > 1) {
    stop("`p` must be between 0.5 and 1.")
  }
  splits <- consensus_tree(trees, p)
  splits <- as.Splits(splits, tipLabels = TipLabels(trees[[1]]))
  summary(splits)
  as.phylo(splits)
}
