test_that("TipsInSplits() family", {
  test <- TipsInSplits(BalancedTree(letters[1:5]))
  expect_identical(test, c("7" = 3L, "8" = 2L, "9" = 2L)[names(test)])
  expect_equal(TipsInSplits(PectinateTree(7), smallest = TRUE),
               c("10" = 2, "11" = 3, "12" = 3, "13" = 2))
  expect_identical(TipsInSplits(PectinateTree(17), keep.names = FALSE), 15:2)

  test <- SplitImbalance(BalancedTree(7))
  expectation <- c("9" = 1L, "10" = 3L, "11" = 3L, "12" = 1L, "13" = 3L)
  expect_identical(test, expectation[names(test)])
})
