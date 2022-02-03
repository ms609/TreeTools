test_that("keep_tip() works", {
  expect_error(keep_tip(BalancedTree(8)$edge[, c(1, 2, 1)], !tabulate(6:4, 8)))
  expect_error(keep_tip(BalancedTree(8)$edge[, 1, drop = FALSE],
                        !tabulate(6:4, 8)))
  expect_error(keep_tip(BalancedTree(8)$edge[, 1], !tabulate(6:4, 8)))
  
  expect_equal(keep_tip(BalancedTree(9)$edge, !tabulate(5:8, 9)),
               matrix(c(6, 7, 8, 9, 9, 8, 7, 6,
                        7, 8, 9, 1, 2, 3, 4, 5), 8, 2))
  
  expect_equal(keep_tip(BalancedTree(8)$edge, !tabulate(6:4, 8)),
               matrix(c(6, 7, 8, 8, 7, 6, 9, 9,
                        7, 8, 1, 2, 3, 9, 4, 5), 8, 2))
  expect_equal(keep_tip(BalancedTree(8)$edge, !tabulate(5:8, 8)),
               BalancedTree(4)$edge)
  
  expect_equal(keep_tip(BalancedTree(8)$edge, !tabulate(3:8, 8)),
               BalancedTree(2)$edge)
  
  expect_equal(keep_tip(unroot(BalancedTree(4))$edge, tabulate(2:4)),
               matrix(c(4, 4, 4, 1, 2, 3), 3, 2))
  
  testTree <- ape::read.tree(text = "(a, ((b, c), ((d, e, f), g)));")
  testEdge <- Preorder(testTree)$edge
})
