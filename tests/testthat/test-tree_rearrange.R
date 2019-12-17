library(ape)

context("Tree rearrangements")

test_that("RootOnNode works", {
  Test <- function (tr, node, rr) {
    if (node <= 8) {
      expect_equal(Preorder(ape::root(tr, outgroup = node, resolve.root = rr)),
                   RootOnNode(tr, node, rr))
    } else {
      expect_equal(Preorder(ape::root(tr, node = node, resolve.root = rr)),
                 RootOnNode(tr, node, rr))
    }
  }

  urt <- UnrootedTreeWithShape(3, 8, letters[1:8])
  Test(urt, 12, TRUE)
  Test(urt, 12, FALSE)
  Test(urt, 4L, TRUE)
  Test(urt, 4L, FALSE)
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
