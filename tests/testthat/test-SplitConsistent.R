test_that("Split contradiction is calculated", {
  pec12 <- as.Splits(PectinateTree(12))
  bal12 <- as.Splits(BalancedTree(12))
  bal8 <- as.Splits(BalancedTree(8))
  expect_error(split_consistent(pec12, list(bal12)),
               "must contain exactly one split")
  expect_error(split_consistent(`attr<-`(pec12[[1]], "nTip", NULL), list(1)),
               "lacks nTip")
  expect_error(split_consistent(`attr<-`(pec12[[1]], "nTip", NULL), list(1)),
               "lacks nTip")
  
  expect_equal(SplitConsistent(pec12[[2]], bal12),
               SplitConsistent(pec12[[2]], list(bal12)))
  expect_error(SplitConsistent(pec12[[2]], "bal12"), "must be a Splits object")
  expect_error(split_consistent(pec12[[2]], bal12), "Not a matrix")
  expect_error(split_consistent(pec12[[2]], list(BalancedTree(12))),
               "not a RawMatrix")
  
  expect_error(split_consistent(bal8[[4]], list(bal12)),
               "different nTip")
  
  expect_equal(split_consistent(pec12[[5]], list(bal12))[[1]],
               rep(TRUE, length(bal12)))
  expect_equal(split_consistent(pec12[[3]], list(bal12))[[1]],
               c(rep(TRUE, 3), FALSE, FALSE, rep(TRUE, 4)))
  expect_equal(SplitConsistent(pec12[[5]], list(bal12, bal12))[[2]],
               rep(TRUE, length(bal12)))
  expect_equal(SplitConsistent(pec12[[5]], list(as.Splits(StarTree(12))))[[1]],
               logical(0))
})
