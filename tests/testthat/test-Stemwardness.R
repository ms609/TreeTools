test_that("Stemwardness functions", {
  bal8 <- BalancedTree(8)
  pec8 <- PectinateTree(8)
  pr3 <- RootTree(pec8, "t3")

  expect_equal(1L, SisterSize(bal8, 3))
  expect_equal(2L, RootNodeDist(bal8, 3))

  expect_equal(5L, SisterSize(pec8, "t3"))
  expect_equal(2L, RootNodeDist(pec8, "t3"))

  expect_equal(7L, SisterSize(pr3, "t3"))
  expect_equal(0L, RootNodeDist(pr3, "t3"))
})
