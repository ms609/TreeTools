context('tree_display')

test_that('ConsensusWithout is robust', {
  expect_equal(BalancedTree(8), ConsensusWithout(BalancedTree(8)))
  expect_equal(BalancedTree(4),
               ConsensusWithout(BalancedTree(8), paste0('t', 5:8)))
  balAndPec <- list(BalancedTree(8), PectinateTree(8))
  t25 <- paste0('t', c(2:5))
  expect_equal(PectinateTree(paste0('t', c('1', 6:8))),
               ConsensusWithout(balAndPec, t25))
  expect_equal(ConsensusWithout(structure(balAndPec, class = 'multiPhylo'), t25),
               ConsensusWithout(balAndPec, t25))
  expect_warning(ConsensusWithout(BalancedTree(8), paste0('t', 1:8)))
})
