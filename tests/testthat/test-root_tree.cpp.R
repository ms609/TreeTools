context("root_tree.cpp")

ApeRoot <- function (tree, root, rr = TRUE) ape::root(tree, root, resolve.root = rr)
Test <- function (tree, root) {
  expect_equal(Preorder(ApeRoot(tree, tree$tip.label[root]))$edge,
               root_on_node(tree$edge, root))
}

test_that('Binary trees are rootable', {
  Test(BalancedTree(9), 3)
  Test(BalancedTree(9), 1)
  Test(PectinateTree(9), 1)
  Test(PectinateTree(9), 7)
})

#test_that('Polytomous trees are rootable', {
  #Test(CollapseNode(BalancedTree(8), 10), 1)
  #Test(CollapseNode(BalancedTree(8), 11), 1)
  #Test(StarTree(8), 1)
#})
