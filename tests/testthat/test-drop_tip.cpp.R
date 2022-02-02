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
  RootingTest <- function (text, tips) {
    tree <- ape::read.tree(text = text)
    expect_equal(
      ape::is.rooted(DropTip(tree, tips)),
      ape::is.rooted(tree)
    )
  }
  
  SeeTree <- function (text, ...) {
    tree <- ape::read.tree(text = text)
    plot(tree)
    nodelabels()
    edgelabels()
    tree$edge
  }
  
  RootingTest("(a, b, (c, d));", 1)
  RootingTest("((a, b), (c, d));", 1)
  RootingTest("(a, b, c, (d, e));", 1)
  RootingTest("(a, (b, (c, d)));", 1)
  
  # Dropping a pair
  RootingTest("((a, b), t2, (t3, t4));", 1:2)
  RootingTest("(a, b, t2, (t3, t4));", 1:2)
  RootingTest("(a, b, (t2, (t3, t4)));", 1:2)
  
  RootingTest("((a, b), (t2, (t3, t4)));", 1:2)
  RootingTest("(a1, a2, b, c, (d, e));", 1:2)

  RootingTest("((a1, a2), (b, (c, d)));", 1:2)
  RootingTest("(a1, a2, (b, (c, d)));", 1:2)
})
