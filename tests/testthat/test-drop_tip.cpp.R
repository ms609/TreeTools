test_that("drop_tip() works", {
  expect_equal(matrix(c(6, 7, 8, 8, 7, 6, 9, 9,
                        7, 8, 1, 2, 3, 9, 4, 5), 8, 2),
               drop_tip(BalancedTree(8)$edge, 6:4))
  expect_equal(BalancedTree(4)$edge,
               drop_tip(BalancedTree(8)$edge, 5:8))

  testTree <- ape::read.tree(text = "(a, ((b, c), ((d, e, f), g)));")
  testEdge <- Preorder(testTree)$edge

  drop_tip(testEdge, 3)
  drop_tip2(testEdge, 3)

  Test <- function (tips) {
    expect_equal(drop_tip(testEdge, tips), drop_tip2(testEdge, tips))
  }
  Test(1)
  Test(2)
  Test(3)
  Test(4)
  Test(5)
  Test(6)
  Test(7)
})
