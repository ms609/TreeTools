context("tree_write.R")

test_that("Write is successful", {
  Test <- function (tree) expect_equal(ape::write.tree(tree), as.Newick(tree))
  Test(BalancedTree(0:7))
  Test(unroot(BalancedTree(0:7)))
  Test(PectinateTree(0:7))
  Test(unroot(PectinateTree(0:7)))
  Test(ape::read.tree(text='(0,1,2,3,4);'))
  Test(ape::read.tree(text='((0,1,2,3,4), ((5, 6), (7, 8, 9)), 10, 11);'))
})
