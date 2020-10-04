#' Distances between each pair of trees
#'
#' @param trees List of trees of class `phylo`.
#' @param Func Function returning a distance between two trees.
#' @param valueLength Integer specifying expected length of the value returned
#' by `Func`.
#' @param \dots Additional arguments to `Func`.
#'
#' @return Matrix detailing distance between each pair of trees.
#' Identical trees are assumed to have zero distance.
#' @examples
#' trees <- list(BalancedTree(8), PectinateTree(8), StarTree(8))
#' TCIDiff <- function (tree1, tree2) {
#'   TotalCopheneticIndex(tree1) - TotalCopheneticIndex(tree2)
#' }
#' PairwiseDistances(trees, TCIDiff, 1)
#' TCIRange <- function (tree1, tree2) {
#'   range(TotalCopheneticIndex(tree1), TotalCopheneticIndex(tree2))
#' }
#' PairwiseDistances(trees, TCIRange, 2)
#' @template MRS
#' @family pairwise tree distances
#' @importFrom stats as.dist
#' @export
PairwiseDistances <- function (trees, Func, valueLength = 1L, ...) {
  ret <- array(0, c(length(trees), length(trees), valueLength))
  for (i in seq_along(trees)) {
    trI <- trees[[i]]
    for (j in i + seq_len(length(trees) - i)) {
      val <- Func(trI, trees[[j]], ...)
      ret[j, i, ] <- unlist(val)
    }
  }

  # Return:
  if (valueLength > 1L) {
    structure(lapply(seq_len(valueLength), function (i) {
      as.dist(ret[, , i], upper = TRUE)
    }), names = names(val))
  } else {
    as.dist(ret[, , 1], upper = TRUE)
  }
}
