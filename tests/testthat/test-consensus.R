test_that("Consensus()", {

  trees <- list(BalancedTree(8), PectinateTree(8))[c(1, 1, 1, 1, 2, 2, 2)]
  trees <- RenumberTips(trees, trees[[1]])
  expect_error(Consensus(trees, 2)) # p too hign
  expect_error(Consensus(trees, 0.2)) # p too low


  expect_equal(RootTree(BalancedTree(8), 't1'),
               RootTree(Consensus(trees, 0.5), 't1'))
  expect_equal(Consensus(trees, 0.999), Consensus(trees, 1))
})
