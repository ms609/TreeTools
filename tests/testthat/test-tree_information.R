context("tree_information.R")

test_that("Phylogenetic info calculated correctly", {
  bal8 <- BalancedTree(8)
  expect_equal(1, TreesMatchingTree(bal8))
  expect_equal(Log2Unrooted(8), PhylogeneticInfo(bal8))
  expect_equal(0L, PhylogeneticInfo(CollapseNode(bal8, 10:15)))

  tr1 <- CollapseNode(bal8, 12:15)
  expect_equal(NUnrooted(4) * NUnrooted(5),
               TreesMatchingTree(tr1))
  expect_equal(Log2Unrooted(8) - log2(45), PhylogeneticInfo(tr1))
})

test_that("Clustering info calculated correctly", {
  bal8 <- BalancedTree(8)
  expect_warning(ClusteringInfo(BalancedTree(8)))
  ciBal8 <- ClusteringInfo(unroot(BalancedTree(8)))
  expect_equal(3L, ciBal8)
  expect_equal(6L, ClusteringInfo(unroot(BalancedTree(16))))
  expect_equal(3L, ClusteringInfo(unroot(PectinateTree(8))))
  expect_gt(ciBal8, ClusteringInfo(CollapseNode(bal8, c(10, 14))))
  expect_gt(ciBal8, ClusteringInfo(CollapseNode(bal8, c(10, 13))))

  expect_equal(2.555, ClusteringInfo(CollapseNode(bal8, c(10, 14))),
               tolerance = 3L)
  expect_equal(2, ClusteringInfo(CollapseNode(bal8, c(10, 13))))
  expect_equal(1.5, ClusteringInfo(CollapseNode(bal8, 10:12)))
  expect_equal(1.58, ClusteringInfo(CollapseNode(bal8, c(10:11, 13))),
               tolerance = 3L)
  expect_equal(2, ClusteringInfo(CollapseNode(PectinateTree(8), c(10, 13))))
  expect_equal(1 + .H(2, 3), ClusteringInfo(CollapseNode(PectinateTree(8), c(10, 14))))

})
