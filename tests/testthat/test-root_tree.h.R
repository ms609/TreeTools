context("root_tree.cpp")

ApeRoot <- function (tree, root, rr = TRUE) ape::root(tree, root, resolve.root = rr)

test_that('Binary trees are rootable', {
  Test <- function (tree, root) {
    expect_equal(Preorder(ApeRoot(tree, tree$tip.label[root]))$edge,
                 root_binary(tree$edge, root))
  }
  Test(BalancedTree(9), 3)
  Test(BalancedTree(9), 1)
  Test(PectinateTree(9), 1)
  Test(PectinateTree(9), 7)
})

test_that('Polytomous trees are rootable', {
  Test <- function (tree, root) {
    expect_equal(Preorder(ApeRoot(tree, tree$tip.label[root])),
                 root_on_node(tree, root))
  }
  bt <- BalancedTree(9)
  pt <- PectinateTree(9)
  Test(CollapseNode(bt, 12), 1)
  Test(CollapseNode(bt, 11), 1)
  Test(CollapseNode(pt, 11), 1)
  Test(CollapseNode(pt, c(11, 12)), 1)
  Test(CollapseNode(pt, c(11, 12)), 5)
  Test(CollapseNode(pt, c(11, 13, 15)), 5)
  Test(StarTree(8), 1)

  t1 <- Preorder(ape::read.tree(text="((10, 7), (6, (8, 11)), (12, (4, (2, 1))), 14, (5, 9, 13), 3);"))
  Test(t1, 1)
  t2 <- Preorder(ape::read.tree(text = "(((2, 4, 5, 7, 9, 10, 12, 13), (1, 14)), (6, (8, 11)), 3);"))
  Test(t2, 1)
})
