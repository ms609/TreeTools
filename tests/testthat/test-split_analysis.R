test_that("TipsInSplits() family", {
  test <- TipsInSplits(BalancedTree(letters[1:5]))
  expect_identical(test, c("7" = 3L, "8" = 2L, "9" = 2L)[names(test)])
  expect_identical(15:2, TipsInSplits(PectinateTree(17), keep.names = FALSE))

  test <- SplitImbalance(BalancedTree(7))
  expectation <- c("9" = 1L, '10' = 3L, '11' = 3L, '12' = 1L, '13' = 3L)
  expect_identical(test, expectation[names(test)])
})
