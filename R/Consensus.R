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
  splits <- as.Splits(consensus_tree(trees, p),
                      tipLabels = TipLabels(trees[[1]]))

  # Return:
  as.phylo(splits)
}


#' @rdname TreeInfo
#'
#' @param trees List of `phylo` objects, optionally with class `multiPhylo`.
#' @param info Abbreviation of 'phylogenetic' or 'clustering', specifying
#' the concept of information to employ.
#' @param check.tips Logical specifying whether to renumber leaves such that
#' leaf numbering is consistent in all trees.
#'
#' @examples
#'
#' # Support-weighted information content of a consensus tree
#' set.seed(0)
#' trees <- list(RandomTree(8), RootTree(BalancedTree(8), 1), PectinateTree(8))
#' cons <- consensus(trees, p = 0.5)
#' p <- SplitFrequency(cons, trees) / length(trees)
#' plot(cons)
#' LabelSplits(cons, signif(SplitwiseInfo(cons, p, sum = FALSE), 4))
#' ConsensusInfo(trees)
#' LabelSplits(cons, signif(ClusteringInfo(cons, p, sum = FALSE), 4))
#' ConsensusInfo(trees, 'clustering')
#' @template MRS
#' @family consensus tree functions
#' @export
ConsensusInfo <- function (trees, info = 'phylogenetic', check.tips = TRUE) {
  mode <- pmatch(tolower(info),
                 c('phylogenetic', 'clustering', 'spic', 'scic')) %% 2
  if (is.na(mode)) {
    stop("`info` must be 'phylogenetic' or 'clustering'")
  }
  if (inherits(trees, 'phylo')) {
    # Convert to length-1 multiphylo object
    trees <- c(trees)
  }
  if (length(trees) == 1L) {
    return((if (mode) SplitwiseInfo else ClusteringInfo)(trees))
  }
  if (check.tips) {
    trees <- RenumberTips(trees, trees[[1]])
  }
  consensus_info(trees, mode == 1L)
}
