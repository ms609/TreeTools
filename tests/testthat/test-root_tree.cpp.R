context("root_tree.cpp")

ApeRoot <- function (tree, root, rr = TRUE) ape::root(tree, root, resolve.root = rr)
Test <- function (tree, root) {
  expect_equal(Preorder(ApeRoot(tree, tree$tip.label[root]))$edge,
               root_on_node(tree$edge, root))
}

test_that('Binary trees are rootable', {
  Test(BalancedTree(8), 3)
  Test(BalancedTree(8), 1)
  Test(PectinateTree(8), 1)
  Test(PectinateTree(8), 7)
})
