context("split_analysis.R")
test_that("TipsInSplits() family", {
  expect_equal(c('8' = 2L, '9' = 2L),
               TipsInSplits(BalancedTree(letters[1:5])))

  expect_equal(c('10' = 3L, '11' = 3L, '12' = 1L, '13' = 3L),
               SplitImbalance(BalancedTree(7)))

})
