test_that('Pairwise distances calculated correctly', {
  nTrees <- 6L
  nTip <- 16L

  set.seed(0)
  trees <- lapply(rep(nTip, nTrees), RandomTree, root = TRUE)
  trees[[1]] <- BalancedTree(nTip)
  trees[[nTrees - 1L]] <- PectinateTree(nTip)
  class(trees) <- 'multiPhylo'

  dists <- PairwiseDistances(trees, phangorn::RF.dist)
  expect_equal(as.integer(phangorn::RF.dist(trees)), as.integer(dists))
})
