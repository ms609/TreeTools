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
  expect_null(ConsensusWithout(BalancedTree(8), paste0('t', 1:8)))

  nasty <- structure(list(edge = structure(
    c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
      5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
    .Dim = c(12, 2)),
    Nnode = 5L,
    tip.label = letters[1:8]),
    class = 'phylo') # Danger: Do not plot!
  expect_equal(Preorder(nasty), ConsensusWithout(nasty))
  expect_equal(DropTip(nasty, 2), ConsensusWithout(nasty, 'b'))

})
