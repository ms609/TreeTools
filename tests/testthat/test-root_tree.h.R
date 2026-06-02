ApeRoot <- function(tree, root, rr = TRUE) ape::root(tree, root, resolve.root = rr)

test_that("Memory leak not encountered", {
  # Example from TreeDist::ClusterTable
  tree1 <- ape::read.tree(text = "(A, (B, (C, (D, E))));");
  tree2 <- ape::read.tree(text = "(A, (B, (D, (C, E))));");
  # as.ClusterTable(tree1) calls: result must be a preordered form of tree1
  expect_equal(Preorder(tree1), root_on_node(tree1, 1))

  # Check for memory leaks...
  root_on_node(RenumberTips(Preorder(tree2), LETTERS[1:5]), 1)[]
  root_on_node(RenumberTips(StarTree(LETTERS[5:1]), LETTERS[1:5]), 1)[]

  expect_error(root_on_node(tree1, 0), "`outgroup` must be a positive integer")
  expect_error(root_on_node(tree1, 999), "`outgroup` exceeds number of nodes")
})

test_that("Big trees don't fail", {
  # 2^14 + 1 is too big for int16; result must be preordered
  expect_equal(root_on_node(PectinateTree(2^14 + 1), 1),
               Preorder(PectinateTree(2^14 + 1)))
})

test_that("root_on_node() always returns preorder (#168)", {
  # Regression: root_on_node() must return PREORDER edges even on its fast path
  # (outgroup already a direct child of the root).  It previously returned the
  # input unchanged there, leaking non-preorder edges into ClusterTable and
  # causing dropped splits / segfaults in consensus_tree() (#168).
  # The earlier already-preorder fixtures don't exercise this: their edge matrix
  # equals Preorder()'s regardless, so they only catch the `order` attribute.
  # We need a tree that is genuinely non-preorder AND has leaf 1 at the root.
  base <- ape::read.tree(text = "(A, (B, (C, (D, E))));") # leaf 1 (A) at root
  nonpre <- base
  nonpre[["edge"]] <- base[["edge"]][c(5, 1, 8, 2, 6, 3, 7, 4), , drop = FALSE]
  attr(nonpre, "order") <- NULL

  # Sanity: input is genuinely non-preorder, but leaf 1 is still a root child
  expect_false(identical(nonpre[["edge"]], Preorder(nonpre)[["edge"]]))
  rootNode <- NTip(nonpre) + 1L
  expect_true(any(nonpre[["edge"]][, 1] == rootNode &
                    nonpre[["edge"]][, 2] == 1L))

  rooted <- root_on_node(nonpre, 1)
  # Fix: edges are reordered to preorder (would equal nonpre's edges pre-fix)
  expect_identical(rooted[["edge"]], Preorder(nonpre)[["edge"]])
  expect_equal(attr(rooted, "order"), "preorder")
  expect_equal(rooted, Preorder(nonpre))
})

test_that("Small trees are rootable", {
  ztt <- ZeroTaxonTree()
  expect_equal(root_on_node(ztt, 1), ztt)
  expect_equal(root_binary(ztt[["edge"]], 1), ztt$edge)
  expect_equal(root_on_node(ztt, 999), ztt)
  
  stt <- SingleTaxonTree()
  expect_equal(root_on_node(stt, 1), stt)
  expect_equal(root_on_node(stt, 2), root_on_node(stt, 1))
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
  expect_error(root_binary(ed9, 9999), "exceeds number of nodes")
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
  expect_error(root_binary(edge, 999), "exceeds number of nodes")
  
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
