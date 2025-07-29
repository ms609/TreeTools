test_that("Split contradiction is calculated", {
  pec12 <- as.Splits(PectinateTree(12))
  bal12 <- as.Splits(BalancedTree(12))
  bal8 <- as.Splits(BalancedTree(8))
  sp_cons <- function(n, h) split_consistent(n, h, FALSE)
  expect_error(sp_cons(pec12, list(bal12)),
               "must contain exactly one split")
  expect_error(sp_cons(pec12[[1]], list(bal8)),
               "Haystack 1 has different nTip than needle")
  expect_error(sp_cons(`attr<-`(pec12[[1]], "nTip", NULL), list(1)),
               "lacks nTip")
  
  expect_equal(SplitConsistent(pec12[[2]], bal12),
               SplitConsistent(pec12[[2]], list(bal12)))
  expect_error(SplitConsistent(pec12[[2]], "bal12"), "must be a Splits object")
  expect_error(SplitConflicts(pec12[[2]], "bal12"), "must be a Splits object")
  expect_error(sp_cons(pec12[[2]], bal12), "Not a matrix")
  expect_error(sp_cons(pec12[[2]], list(BalancedTree(12))),
               "not a RawMatrix")
  
  expect_error(sp_cons(bal8[[4]], list(bal12)),
               "different nTip")
  
  expect_equal(sp_cons(pec12[[5]], list(bal12))[[1]],
               rep(TRUE, length(bal12)))
  expect_equal(sp_cons(pec12[[3]], list(bal12))[[1]],
               c(rep(TRUE, 3), FALSE, FALSE, rep(TRUE, 4)))
  expect_equal(split_consistent(pec12[[3]], list(bal12), TRUE)[[1]],
               !c(rep(TRUE, 3), FALSE, FALSE, rep(TRUE, 4)))
  expect_equal(SplitConsistent(pec12[[5]], list(bal12, bal12))[[2]],
               rep(TRUE, length(bal12)))
  expect_equal(SplitConflicts(pec12[[5]], as.Splits(StarTree(12)))[[1]],
               logical(0))
  expect_equal(SplitConflicts(pec12[[2]], list(bal12, pec12)),
               lapply(SplitConsistent(pec12[[2]], list(bal12, pec12)), `!`))
})
