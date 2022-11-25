test_that("as.matrix(phylo)", {
  bal8 <- BalancedTree(8)
  sp8 <- as.Splits(bal8)
  expect_equal(as.matrix(bal8), 1 * as.logical(sp8))
  expect_equal(as.logical(bal8), as.logical(sp8))
  expect_equal(as.matrix(sp8), as.logical(sp8))
}
