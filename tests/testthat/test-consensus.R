test_that("Consensus()", {

  ApeTest <- function (tr, p = 1) {
    # plot(Consensus(tr))
    # plot(consensus(tr))
    expect_true(all.equal(RootTree(Consensus(tr, p = p), 1),
                          RootTree(ape::consensus(tr, p = p), 1)))
  }

  trees <- list(BalancedTree(8), PectinateTree(8))[c(1, 1, 1, 1, 2, 2, 2)]
  trees <- RenumberTips(trees, trees[[1]])
  expect_error(Consensus(trees, 2)) # p too hign
  expect_error(Consensus(trees, 0.2)) # p too low

  ApeTest(trees)
  ApeTest(trees, 0.5)
  expect_equal(Consensus(trees, 0.999), Consensus(trees, 1))

  ApeTest(as.phylo(0:2, 8))
  ApeTest(as.phylo(0:250, 8))
  ApeTest(as.phylo(0:250, 80))

})
