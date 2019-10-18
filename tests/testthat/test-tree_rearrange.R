library(ape)

context("Tree rearrangements")

test_that("CollapseNodes works", {
  tree8  <- read.tree(text="(((a, (b, (c, d))), (e, f)), (g, h));")
  expect_error(CollapseNode(1:5, tree8))
  expect_error(CollapseNode(tree8, 1))
  
  suppressWarnings(RNGversion("3.5.0")) # Until we can require R3.6.0
  set.seed(1)
  
  tree <- rtree(7)
  expect_equal(tree, CollapseNode(tree, integer(0)))
  
  no1213 <- CollapseNode(tree, c(12, 13))
  expect_equal(no1213$edge, matrix(c(8, 9, 9, 8, 10, 11, 11, 10, 10, 10, 
                                      9, 1, 2, 10, 11, 3:7), ncol=2))
  el <- tree$edge.length
  expect_equal(no1213$edge.length, c(el[1:7], el[8] + c(c(el[9] + el[10:11]), el[12])))
  
  no11 <- CollapseEdge(tree, 5L)
  expect_equal(no11$edge, matrix(c(8, 9, 9, 8, 10, 10, 10, 11, 12, 12, 11,
                                   9, 1, 2, 10, 3, 4, 11, 12, 5:7), ncol=2))

})
