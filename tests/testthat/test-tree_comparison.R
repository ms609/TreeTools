test_that("Pairwise distances calculated correctly", {
  nTrees <- 6L
  nTip <- 16L

  set.seed(0)
  trees <- lapply(rep(nTip, nTrees), RandomTree, root = TRUE)
  trees[[1]] <- BalancedTree(nTip)
  trees[[nTrees - 1L]] <- PectinateTree(nTip)
  class(trees) <- "multiPhylo"

  # From example
  TCIRange <- function(tree1, tree2) {
   range(TotalCopheneticIndex(tree1), TotalCopheneticIndex(tree2))
  }
  tciPairs <- PairwiseDistances(trees, TCIRange, 2)
  expect_equal(length(tciPairs), 2)
  expect_equal(as.matrix(tciPairs[[1]])[3, 6],
               TCIRange(trees[[3]], trees[[6]])[1])
  
  skip_if_not_installed("phangorn")
  trees <- reorder(trees, "cladewise")
  dists <- PairwiseDistances(trees, phangorn::RF.dist)
  expect_equal(as.integer(phangorn::RF.dist(trees)), as.integer(dists))
})
