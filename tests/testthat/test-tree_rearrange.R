library(ape)

context("Tree rearrangements")

test_that("RootOnNode works", {
  TestTip <- function (tr, node, rr) {
    expect_equal(Preorder(ape::root(tr, outgroup = node, resolve.root = rr)),
                 RootOnNode(tr, node, rr))
  }
  TestInternal <- function (tr, node, rr) {
    expect_equal(Preorder(ape::root(tr, node = node, resolve.root = rr)),
                 RootOnNode(tr, node, rr))
  }

  TestInternal(PectinateTree(8), 14, TRUE)
  TestInternal(PectinateTree(8), 14, FALSE)
  TestTip(PectinateTree(8), 4L, TRUE)
  TestTip(PectinateTree(8), 4L, FALSE)

  TestInternal(BalancedTree(8L), 11L, TRUE)
  TestInternal(BalancedTree(8L), 11L, FALSE)

  urt <- UnrootedTreeWithShape(3, 8, letters[1:8])
  TestInternal(urt, 12, TRUE)
  TestInternal(urt, 12, FALSE)
  TestTip(urt, 4L, TRUE)
  TestTip(urt, 4L, FALSE)

  # Children of an unresolved root
  TestInternal(urt, 10, TRUE)
  TestTip(urt, 1, TRUE)

  expect_equal(PectinateTree(8), RootOnNode(PectinateTree(8), 9L, TRUE))
  expect_equal(unroot(PectinateTree(8)), RootOnNode(PectinateTree(8), 9L, FALSE))
  expect_equal(urt, RootOnNode(urt, 9L, FALSE))
  expect_equal(Preorder(EnforceOutgroup(urt, letters[1:2])),
               RootOnNode(urt, 9L, TRUE))

})

test_that("CollapseNodes works", {
  tree8  <- read.tree(text="(((a, (b, (c, d))), (e, f)), (g, h));")
  expect_error(CollapseNode(1:5, tree8))
  expect_error(CollapseNode(tree8, 1))

  tree <- as.phylo(123, 7)
  tree$edge.length <- 12:1
  expect_equal(tree, CollapseNode(tree, integer(0)))

  exp1213 <- matrix(c(8, 8, 9, 10, 10, 9, 11, 11, 11, 11,
                      1, 9, 10, 2, 4, 11, 3, 7, 6, 5), ncol=2)
  no1213 <- CollapseNode(tree, c(12, 13))
  expect_equal(exp1213, no1213$edge)

  el <- tree$edge.length
  expect_equal(no1213$edge.length, c(el[1:6],
                                     el[7] + c(c(el[8] + el[9:10]), el[11]),
                                     el[12]))

  no11 <- CollapseEdge(tree, c(7, 8))
  expect_equal(exp1213, no11$edge)

})
