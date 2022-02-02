test_that("drop_tip() works", {
  expect_error(drop_tip(BalancedTree(8)$edge[, c(1, 2, 1)], 6:4))
  expect_error(drop_tip(BalancedTree(8)$edge[, 1, drop = FALSE], 6:4))
  expect_error(drop_tip(BalancedTree(8)$edge[, 1], 6:4))

  expect_equal(matrix(c(6, 7, 8, 8, 7, 6, 9, 9,
                        7, 8, 1, 2, 3, 9, 4, 5), 8, 2),
               drop_tip(BalancedTree(8)$edge, 6:4))
  expect_equal(BalancedTree(4)$edge,
               drop_tip(BalancedTree(8)$edge, 5:8))

  testTree <- ape::read.tree(text = "(a, ((b, c), ((d, e, f), g)));")
  testEdge <- Preorder(testTree)$edge
})

test_that("drop_tip() retains rootedness", {
  Edge <- function (text) ape::read.tree(text = text)$edge
  
  expect_equal(drop_tip(Edge("(a, b, (c, d));"), 1),
               Edge("(a, b, c);"))
  
  expect_equal(drop_tip(Edge("(a, b, c, (d, e));"), 1),
               Edge("(a, b, (c, d));"))
  
  expect_equal(drop_tip(Edge("(a, (b, (c, d)));"), 1),
               Edge("(b, (c, d));"))
  
  
  # Dropping a pair
  
  expect_equal(drop_tip(Edge("((a, b), t2, (t3, t4));"), 1:2),
               Edge("(a, b, c);"))
  
  expect_equal(drop_tip(Edge("(a, b, t2, (t3, t4));"), 1:2),
               Edge("(a, b, c);"))
  
  expect_equal(drop_tip(Edge("(a, b, (t2, (t3, t4)));"), 1:2),
               Edge("(a, b, c);"))
  
  expect_equal(drop_tip(Edge("((a, b), (t2, (t3, t4)));"), 1:2),
               Edge("(a, b, c);"))
  
  
  expect_equal(drop_tip(Edge("(a1, a2, b, c, (d, e));"), 1:2),
               Edge("(a, b, (c, d));"))
  
  expect_equal(drop_tip(Edge("((a1, a2), (b, (c, d)));"), 1:2),
               Edge("(b, (c, d));"))
  
  expect_equal(drop_tip(Edge("(a1, a2, (b, (c, d)));"), 1:2),
               Edge("(b, (c, d));"))
  
})
