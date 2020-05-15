context("tree_write.R")

test_that("Write is successful", {
  Test <- function (tree) expect_equal(ape::write.tree(tree), as.Newick(tree))
  Test(BalancedTree(0:7))
  Test(unroot(BalancedTree(0:7)))
  Test(PectinateTree(0:7))
  Test(unroot(PectinateTree(0:7)))
  Test(ape::read.tree(text='(0,1,2,3,4);'))
  Test(ape::read.tree(text='((0,1,2,3,4), ((5, 6), (7, 8, 9)), 10, 11);'))

  nasty <- structure(list(edge = structure(
    c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
      5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
    .Dim = c(12, 2)),
    Nnode = 5L,
    tip.label = 0:7),
    class = 'phylo') # Danger: Do not plot!
  Test(nasty)
})
