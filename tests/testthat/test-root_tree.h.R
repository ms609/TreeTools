ApeRoot <- function(tree, root, rr = TRUE) ape::root(tree, root, resolve.root = rr)

test_that("Memory leak not encountered", {
  # Example from TreeDist::ClusterTable
  tree1 <- ape::read.tree(text = "(A, (B, (C, (D, E))));");
  tree2 <- ape::read.tree(text = "(A, (B, (D, (C, E))));");
  # as.ClusterTable(tree1) calls:
  expect_equal(tree1, root_on_node(tree1, 1))

  # Check for memory leaks...
  root_on_node(RenumberTips(Preorder(tree2), LETTERS[1:5]), 1)[]
  root_on_node(RenumberTips(StarTree(LETTERS[5:1]), LETTERS[1:5]), 1)[]

  expect_error(root_on_node(tree1, 0), "`outgroup` must be a positive integer")
  expect_error(root_on_node(tree1, 999), "`outgroup` exceeds number of nodes")
})

test_that("Big trees don't fail", {
  # 2^14 + 1 is too big for int16
  expect_equal(root_on_node(PectinateTree(2^14 + 1), 1),
               PectinateTree(2^14 + 1))
})

test_that("Binary trees are rootable", {
  Test <- function(tree, root) {
    expect_equal(Preorder(ApeRoot(tree, tree$tip.label[root]))$edge,
                 root_binary(tree$edge, root))
  }
  Test(BalancedTree(9), 3)
  Test(BalancedTree(9), 1)
  Test(PectinateTree(9), 1)
  Test(PectinateTree(9), 7)
  ed9 <- PectinateTree(9)$edge
  expect_equal(root_binary(ed9, 10), ed9)
  expect_equal(root_binary(ed9, 1), ed9)
})

test_that("Polytomous trees are rootable", {
  Test <- function(tree, root) {
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
  t1 <- Preorder(ape::read.tree(
    text = "((10, 7), (6, (8, 11)), (12, (4, (2, 1))), 14, (5, 9, 13), 3);"))
  Test(t1, 1)
  t2 <- Preorder(ape::read.tree(
    text = "(((2, 4, 5, 7, 9, 10, 12, 13), (1, 14)), (6, (8, 11)), 3);"))
  Test(t2, 1)
})

test_that("Rooted trees report preorder accurately", {
  set.seed(1)
  nTips <- 8
  edge <- do.call(cbind,
                  RenumberEdges(.RandomParent(nTips),
                                seq_len(nTips + nTips - 2L)))
  
  expect_preorder <- function(x) {
    expect_equal(x, Preorder(`attr<-`(x, "order", "unknown")))
  }
  
  # Check that we are in preorder
  expect_preorder(edge)
  
  expect_preorder(root_binary(edge, 2))
  expect_preorder(root_binary(edge, 6))
  
  tree <- structure(list(edge = edge,
                         Nnode = nTips - 1L,
                         tip.label = TipLabels(nTips)),
                    order = "preorder",
                    class = "phylo")
  
  rootNode <- nTips + 1L
  expect_preorder(root_on_node(tree, rootNode))
  deepNode <- 2 * nTips - 2
  expect_preorder(root_on_node(tree, deepNode))
  
  weighted <- tree
  weights <- seq_len(dim(edge)[[1]])
  weighted[["edge.length"]] <- weights
  expect_equal(sort(root_on_node(weighted, rootNode)[["edge.length"]]),
                    weighted[["edge.length"]])
  expect_preorder(root_on_node(weighted, rootNode))
  expect_preorder(root_on_node(weighted, deepNode))
})
