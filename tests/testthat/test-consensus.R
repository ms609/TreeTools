test_that("Consensus()", {

  ApeTest <- function (tr) {
    # plot(Consensus(tr))
    # plot(consensus(tr))
    expect_true(all.equal(Consensus(tr), ape::consensus(tr)))
  }

  trees <- list(BalancedTree(8), PectinateTree(8))[c(1, 1, 1, 1, 2, 2, 2)]
  trees <- RenumberTips(trees, trees[[1]])
  expect_error(Consensus(trees, 2)) # p too hign
  expect_error(Consensus(trees, 0.2)) # p too low

  ApeTest(trees)
  expect_equal(RootTree(BalancedTree(8), 't1'),
               RootTree(Consensus(trees, 0.5), 't1'))
  expect_equal(Consensus(trees, 0.999), Consensus(trees, 1))
  ApeTest(as.phylo(0:2, 8))
  par(mfrow=2:1, mar = rep(0, 4))
  plot(Consensus(as.phylo(0:2, 8))); plot(consensus(as.phylo(0:2, 8))); legend('bottomleft', 'APE')
  ApeTest(as.phylo(0:250, 8))
  ApeTest(as.phylo(0:250, 80))
  ApeTest(as.phylo(0:250, 800))
})
