test_that("EdgeAncestry() works", {
  tree <- BalancedTree(10)
  edge <- tree$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  expect_equal(which(EdgeAncestry(4L, parent, child)), 1:3)
  expect_equal(which(EdgeAncestry(2, parent, child)), 1)
  expect_equal(which(EdgeAncestry(1, parent, child)), integer(0))
  expect_equal(EdgeAncestry(10, parent, child), logical(18))
})

test_that("MRCA() works", {
  bal7 <- BalancedTree(7)
  allAnc <- AllAncestors(bal7$edge[, 1], bal7$edge[, 2])
  expect_equal(9, MRCA(1, 4, allAnc))
  expect_equal(8, MRCA(1, 6, allAnc))
  expect_equal(8, MRCA(1, 7, allAnc))
  expect_equal(1, MRCA(1, 1, allAnc))
  expect_equal(9, MRCA(1, 11, allAnc))
})
