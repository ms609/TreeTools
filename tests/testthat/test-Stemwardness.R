context("Stemwardness.R")

test_that("Stemwardness()", {
  bal8 <- BalancedTree(8)
  pec8 <- PectinateTree(8)
  expect_equal(c(sisterSize = 1L, rootNodeDist = 2L), Stemwardness(bal8, 3))
  expect_equal(c(sisterSize = 5L, rootNodeDist = 2L), Stemwardness(pec8, 't3'))
  expect_equal(c(sisterSize = 7L, rootNodeDist = 0L),
               Stemwardness(RootTree(pec8, 't3'), 't3'))
})
