ApeRoot <- function (tree, root, rr = TRUE) ape::root(tree, root, resolve.root = rr)

test_that("Memory leak not encountered", {
  # Example from TreeDist::ClusterTable
  tree1 <- ape::read.tree(text = "(A, (B, (C, (D, E))));");
  tree2 <- ape::read.tree(text = "(A, (B, (D, (C, E))));");
  # as.ClusterTable(tree1) calls:
  expect_equal(tree1, root_on_node(tree1, 1))

  # Check for memory leaks...
  root_on_node(RenumberTips(Preorder(tree2), LETTERS[1:5]), 1)[]
  root_on_node(RenumberTips(StarTree(LETTERS[5:1]), LETTERS[1:5]), 1)[]
})

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
  Test(CollapseNode(bt, 12), 3)
  Test(CollapseNode(bt, 11), 1)
  Test(CollapseNode(pt, 11), 1)
  Test(CollapseNode(pt, c(11, 12)), 1)
  Test(CollapseNode(pt, c(11, 12)), 3)
  Test(CollapseNode(pt, c(11, 12)), 5)
  Test(CollapseNode(pt, c(11, 13, 15)), 5)
  Test(CollapseNode(pt, c(11:13, 15)), 9)
  Test(StarTree(8), 1)

  # Day 1985 examples
  t1 <- Preorder(ape::read.tree(text="((10, 7), (6, (8, 11)), (12, (4, (2, 1))), 14, (5, 9, 13), 3);"))
  Test(t1, 1)
  t2 <- Preorder(ape::read.tree(text = "(((2, 4, 5, 7, 9, 10, 12, 13), (1, 14)), (6, (8, 11)), 3);"))
  Test(t2, 1)
})
