context("TotalCopheneticIndex.R")

test_that("Trees from Mir et al. 2013 are scored correctly", {
  Tree <- function (text) ape::read.tree(text = text)
  expect_equal(0L,  TotalCopheneticIndex(Tree('(1,2,3,4,5);')))
  expect_equal(1L,  TotalCopheneticIndex(Tree('((1,2),3,4,5);')))
  expect_equal(2L,  TotalCopheneticIndex(Tree('((1,2),(3,4),5);')))
  expect_equal(3L,  TotalCopheneticIndex(Tree('((1,2,3),4,5);')))
  expect_equal(4L,  TotalCopheneticIndex(Tree('(((1,2),3),4,5);')))
  expect_equal(4L,  TotalCopheneticIndex(Tree('((1,2,3),(4,5));')))
  expect_equal(5L,  TotalCopheneticIndex(Tree('(((1,2),3),(4,5));')))
  expect_equal(6L,  TotalCopheneticIndex(Tree('((1,2,3,4),5);')))
  expect_equal(7L,  TotalCopheneticIndex(Tree('(((1,2),3,4),5);')))
  expect_equal(8L,  TotalCopheneticIndex(Tree('(((1,2),(3,4)),5);')))
  expect_equal(9L,  TotalCopheneticIndex(Tree('(((1,2,3),4),5);')))
  expect_equal(10L, TotalCopheneticIndex(Tree('((((1,2),3),4),5);')))
  expect_equal(10L, TotalCopheneticIndex(PectinateTree(5)))
  expect_equal(10L, TotalCopheneticIndex(PectinateTree(5)))


  nasty <- structure(list(edge = structure(
    c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
      5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
    .Dim = c(12, 2)),
    Nnode = 5L,
    tip.label = letters[1:8]),
    class = 'phylo')
  expect_equal(28L, TotalCopheneticIndex(nasty))

  expect_equal(TCIContext(BalancedTree(5)), TCIContext(5L))
})
