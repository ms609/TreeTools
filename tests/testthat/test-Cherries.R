test_that("Cherries works", {
  expect_equal(Cherries(BalancedTree(8)), 4L)
  expect_equal(Cherries(PectinateTree(8)), 1L)
  expect_equal(Cherries(UnrootTree(PectinateTree(8))$edge, 8L), 2L)
  expect_error(Cherries(matrix(4, 4, 4)), "edge matrix")
  expect_error(Cherries(1:3), "edge matrix")
  
  expect_error(n_cherries_wrapper(1:2, 1:3, 4), "same length")
  
  expect_no_error(Cherries(BalancedTree(144)))
})
