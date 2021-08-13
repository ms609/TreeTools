test_that("Consensus() errors", {
  bal8 <- BalancedTree(8)
  oneLeaf <- Consensus(list(DropTip(bal8, 1:7))[c(1, 1, 1)])
  expect_equal(class(oneLeaf), 'phylo')
  expect_equal(oneLeaf$Nnode, 0)
  expect_identical(Consensus(list(DropTip(bal8, 1:8))[c(1, 1, 1)]),
                   DropTip(bal8, 1:8))
  expect_equal(expect_warning(Consensus(list(PectinateTree(6), PectinateTree(8)))),
               PectinateTree(6))
})

test_that("Consensus()", {

  ApeTest <- function (tr, p = 1) {
    # plot(Consensus(tr))
    # plot(consensus(tr))
    if (!expect_true(all.equal(RootTree(Consensus(tr, p = p), 1),
                               RootTree(ape::consensus(tr, p = p), 1)))) {
      dput(RootTree(Consensus(tr, p = p), 1)$edge)
      dput(RootTree(ape::consensus(tr, p = p), 1)$edge)
    }
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

test_that("ConsensusWithout()", {
  tr <- as.phylo(0:6, 16)
  expect_identical(ConsensusWithout(tr, paste0('t', 1:16)),
                   DropTip(tr[[1]], 1:16))
  expect_identical(ConsensusWithout(tr, paste0('t', 1:4)),
                   Consensus(DropTip(tr, 1:4)))
})
