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
  expect_warning(ClusterInfo(BalancedTree(8)))
  ciBal8 <- ClusterInfo(unroot(BalancedTree(8)))
  expect_equal(3L, ciBal8)
  expect_equal(6L, ClusterInfo(unroot(BalancedTree(16))))
  expect_equal(3L, ClusterInfo(unroot(PectinateTree(8))))
  halfNHalf <- ape::drop.tip(ape::bind.tree(BalancedTree(8), PectinateTree(9)), 9)

  ClusterInfo(halfNHalf)

  expect_gt(ciBal8, ClusterInfo(CollapseNode(bal8, c(10, 14))))
  expect_gt(ciBal8, ClusterInfo(CollapseNode(bal8, c(10, 13))))

  expect_equal(2.555, ClusterInfo(CollapseNode(bal8, c(10, 14))),
               tolerance = 3L)
  expect_equal(2, ClusterInfo(CollapseNode(bal8, c(10, 13))))
  expect_equal(1.5, ClusterInfo(CollapseNode(bal8, 10:12)))
  expect_equal(1.58, ClusterInfo(CollapseNode(bal8, c(10:11, 13))),
               tolerance = 3L)
  expect_equal(2, ClusterInfo(CollapseNode(PectinateTree(8), c(10, 13))))
  expect_equal(1 + .H(2, 3), ClusterInfo(CollapseNode(PectinateTree(8), c(10, 14))))

})
