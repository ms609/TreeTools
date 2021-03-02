context('tree_display')

test_that('ConsensusWithout() is robust', {
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

test_that("SortTree() works", {
  expect_equal(matrix(c(7:10, 10, 9, 11, 11, 8, 7:10, 3:4, 11, 1:2, 5:6), 10),
               SortTree(as.phylo(10, 6))$edge)
  expect_error(#TODO sort unrooted trees,
               SortTree(UnrootTree(PectinateTree(5)))$edge)
  expect_equal('cladewise', attr(SortTree(PectinateTree(5)), 'order'))
})

test_that("SortTree.multiPhylo()", {
  t1 <- as.phylo(123, 12)
  t2 <- as.phylo(921, 12)

  expect_identical(structure(list(SortTree(t1), SortTree(t2)),
                             class = 'multiPhylo'),
                   SortTree(c(t1, t2)))
})
