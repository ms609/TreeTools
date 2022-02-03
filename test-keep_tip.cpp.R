test_that("keep_tip() works", {
  expect_error(keep_tip(BalancedTree(8)$edge[, c(1, 2, 1)], !tabulate(6:4, 8)))
  expect_error(keep_tip(BalancedTree(8)$edge[, 1, drop = FALSE],
                        !tabulate(6:4, 8)))
  expect_error(keep_tip(BalancedTree(8)$edge[, 1], !tabulate(6:4, 8)))
  
  expect_equal(matrix(c(6, 7, 8, 8, 7, 6, 9, 9,
                        7, 8, 1, 2, 3, 9, 4, 5), 8, 2),
               keep_tip(BalancedTree(9)$edge, !tabulate(5:8, 9)))
  
  expect_equal(matrix(c(6, 7, 8, 8, 7, 6, 9, 9,
                        7, 8, 1, 2, 3, 9, 4, 5), 8, 2),
               keep_tip(BalancedTree(8)$edge, !tabulate(6:4, 8)))
  expect_equal(BalancedTree(4)$edge,
               keep_tip(BalancedTree(8)$edge, !tabulate(5:8, 8)))
  
  testTree <- ape::read.tree(text = "(a, ((b, c), ((d, e, f), g)));")
  testEdge <- Preorder(testTree)$edge
})
