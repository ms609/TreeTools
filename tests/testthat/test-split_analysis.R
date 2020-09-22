context("split_analysis.R")
test_that("TipsInSplits() family", {
  expect_identical(c('8' = 2L, '9' = 2L),
                   TipsInSplits(BalancedTree(letters[1:5])))
  expect_identical(15:2, TipsInSplits(PectinateTree(17), keep.names = FALSE))

  expect_identical(c('10' = 3L, '11' = 3L, '12' = 1L, '13' = 3L),
                   SplitImbalance(BalancedTree(7)))

})
