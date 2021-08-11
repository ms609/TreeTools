test_that("Consensus()", {

  trees <- list(BalancedTree(8), PectinateTree(8))[c(1, 1, 1, 1, 2, 2, 2)]
  trees <- RenumberTips(trees, trees[[1]])
  expect_error(Consensus(trees, 2)) # p too hign
  expect_error(Consensus(trees, 0.2)) # p too low


  expect_equal(RootTree(BalancedTree(8), 't1'),
               RootTree(Consensus(trees, 0.5), 't1'))
  expect_equal(Consensus(trees, 0.999), Consensus(trees, 1))
})

test_that("ConsensusInfo() is robust", {
  trees <- list(read.tree(text = '(a, (b, (c, (d, (e, X)))));'),
                read.tree(text = '((a, X), (b, (c, (d, e))));'))
  expect_equal(0, ConsensusInfo(trees, 'cl'))
  expect_error(ConsensusInfo(trees, 'ERROR'))
  # multiPhylo vs phylo
  expect_equal(ConsensusInfo(trees[1]), ConsensusInfo(trees[[1]]))
})

test_that("ConsensusInfo() generates correct value", {
  trees <- list(read.tree(text = "((a, b), (c, d));"),
                read.tree(text = "((a, c), (b, d));"),
                read.tree(text = "((a, d), (c, b));"))
  expect_equal(0, ConsensusInfo(trees))
  expect_equal(0, ConsensusInfo(trees, 'cl'))
  expect_equal(log2(3), ConsensusInfo(trees[1]))
  expect_equal(4, ConsensusInfo(trees[1], 'cl'))
  expect_equal(log2(3), ConsensusInfo(trees[c(1, 1)]))
  expect_equal(4, ConsensusInfo(trees[c(1, 1)], 'cl'))

  Entropy <- function (...) {
    p <- c(...)
    p <- p[p > 0]
    -sum(p * log2(p))
  }

  expect_equal(Entropy(c(1, 1, 1) / 3) - Entropy(c(1/2, 1/2, 9)/10),
               ConsensusInfo(trees[c(rep(1, 9), 2)]))
})
