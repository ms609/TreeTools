context("tree_information.R")

test_that("Cladistic info calculated correctly", {
  bal8 <- BalancedTree(8)
  expect_equal(1, TreesMatchingTree(bal8))
  expect_equal(0, LnTreesMatchingTree(bal8))
  expect_equal(Log2Unrooted(8), CladisticInfo(bal8))
  expect_equal(0L, CladisticInfo(StarTree(8)))
  expect_equal(rep(CladisticInfo(PectinateTree(6)), 4),
               CladisticInfo(as.phylo(1:4, 6)))

  tr1 <- CollapseNode(bal8, 12:15)
  expect_equal(NUnrooted(4) * NUnrooted(5),
               TreesMatchingTree(tr1))
  expect_equal(Log2Unrooted(8) - log2(45), CladisticInfo(tr1))
})
