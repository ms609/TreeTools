context("tree_properties.R")

test_that('EdgeAncestry works', {
  tree <- BalancedTree(10)
  edge <- tree$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  expect_equal(1:3, which(EdgeAncestry(4L, parent, child)))
  expect_equal(1, which(EdgeAncestry(2, parent, child)))
  expect_equal(integer(0), which(EdgeAncestry(1, parent, child)))
  expect_equal(logical(18), EdgeAncestry(10, parent, child))
})

test_that('Root node can be found', {
  rooted <- BalancedTree(8)
  postorder <- Postorder(rooted)
  unrooted <- UnrootTree(rooted)
  expect_true(TreeIsRooted(rooted))
  expect_false(TreeIsRooted(unrooted))
  expect_equal(9L, RootNode(rooted$edge))
  expect_equal(9L, RootNode(rooted))
  expect_equal(9L, RootNode(as.phylo(1337L, 8L)))
  expect_equal(rep(9L, 3L),
               RootNode(c(unrooted, Preorder(postorder), postorder)))
  expect_equal(rep(9L, 3L),
               RootNode(list(Cladewise(postorder),
                             ApePostorder(rooted),
                             Pruningwise(postorder))))
  expect_warning(RootNode(matrix(1:4, 2)))
})

test_that('NodeOrder() works', {
  expect_equivalent(c(2L, rep(3L, 6)),
                    as.integer(NodeOrder(BalancedTree(8), internalOnly = TRUE)))
  expect_equivalent(c(rep(1L, 8), 2L, rep(3L, 6)),
                    as.integer(NodeOrder(BalancedTree(8), internalOnly = FALSE)))
  expect_equivalent(rep(2L, 7), as.integer(NodeOrder(BalancedTree(8), FALSE)))
  tree <- CollapseNode(BalancedTree(8), 12:15)
  expect_equivalent(c(`9` = 5L, `10` = 4L, `11` = 3L),
                    as.integer(NodeOrder(tree$edge, internalOnly = TRUE)))
  expect_equal(list(NodeOrder(tree), NodeOrder(tree)), NodeOrder(c(tree, tree)))
  expect_equal(NodeOrder(c(tree, tree)), NodeOrder(list(tree, tree)))
})

test_that('Rooting and partition counting', {
  set.seed(0)

  expect_false(TreeIsRooted(ape::read.tree(text='(a, b, c);')))
  expect_true(TreeIsRooted(ape::read.tree(text='(a, (b, c));')))

  tree8 <- RandomTree(8L, TRUE)
  expect_true(TreeIsRooted(tree8))
  expect_false(TreeIsRooted(UnrootTree(tree8)))

  expect_equal(5L, NPartitions(8))
  expect_equal(5L, NPartitions(tree8))
  expect_equal(5L, NPartitions(UnrootTree(tree8)))

  expect_equal(0L, NPartitions(ape::read.tree(text='(A, ((B, C)));')))

  expect_equal(c(5L, 5L), NPartitions(list(tree8, UnrootTree(tree8))))
  expect_equal(c(5L, 5L), NPartitions(c(8, 8)))
  expect_error(NPartitions('not a tree'))
})

test_that("NTip() works", {
  Test <- function (n) {
    tr <- BalancedTree(n)
    pec <- PectinateTree(n)
    expect_identical(n, NTip(tr))
    expect_identical(n, NTip(tr$edge))
    expect_identical(n, NTip(Postorder(tr$edge)))
    expect_identical(n, NTip(list(tr)))
    expect_identical(rep(n, 2L), NTip(list(tr, tr)))
    expect_identical(rep(n, 3L),
                     NTip(structure(list(tr, tr, pec), class = 'multiPhylo')))
    expect_error(NTip(matrix('', n, 2)))
  }
  Test(1L)
  Test(8L)
  Test(64L)
  Test(2000L)

})

test_that('Edge distances are calculated correctly', {
  tree <- ape::read.tree(text = '(((a, b), (c, d)), (e, (f, (g, (h, i)))));')
  ed <- EdgeDistances(tree)

  expect_equal(ed, t(ed)) # Symmetry
  expect_equal(c(4, 5, 6, 6,
                 5, 6, 6, 4,
                 4, 3, 3, 2,
                 2, 1, 1, 0),
               ed[, 16])
  expect_equal(c(1, 2, 3, 3, 2, 3, 3, 1, 1, 0, 1, 1, 2, 2, 3, 3), ed[, 10])
  expect_equal(c(1, 0, 1, 1, 1, 2, 2, 1, 2, 2, 3, 3, 4, 4, 5, 5), ed[, 2])

  tree <- ape::read.tree(text='(a, (b, (c, (d, (((e, f), g), (h, (i, j)))))));')
  ed <- EdgeDistances(tree)
  expect_equal(ed, t(ed)) # Symmetry
  expect_equal(c(6, 6, 6, 5, 5, 4, 4, 3, 2, 1, 0, 1, 2, 3, 4, 4, 5, 5),
               ed[11, ])
})

test_that("Node depths calculated correctly", {
  #par(mar=rep(0.4, 4))
  expect_equal(c(rep(0, 20), 19:1), NodeDepth(PectinateTree(20)))
  expect_equal(c(rep(0, 20), rep(1, 19)),
               NodeDepth(PectinateTree(20), shortest = TRUE))
  expect_equal(c(rep(0, 20), 1:9, 9:1), NodeDepth(UnrootTree(PectinateTree(20))))
  expect_equal(c(rep(0, 20), rep(1, 18)),
               NodeDepth(UnrootTree(PectinateTree(20)), shortest = TRUE))

  tree <- BalancedTree(20)
  expect_equal(c(5,4,3,2,1,1,3,2,1,1,4,3,2,1,1,3,2,1,1),
               NodeDepth(tree, includeTips = FALSE))
  expect_equal(c(4,3,2,1,1,1,2,1,1,1,3,2,1,1,1,2,1,1,1),
               NodeDepth(tree, shortest = TRUE, includeTips = FALSE))

  tree <- UnrootTree(tree)
  expect_equal(c(4,3,2,1,1,3,2,1,1,4,3,2,1,1,3,2,1,1),
               NodeDepth(tree, includeTips = FALSE))
  expect_equal(c(3,2,1,1,1,2,1,1,1,3,2,1,1,1,2,1,1,1),
               NodeDepth(tree, shortest = TRUE, includeTips = FALSE))

  tree <- UnrootTree(RootTree(BalancedTree(20), 't10'))
  #plot(tree)
  #nodelabels(NodeDepth(tree, F, FALSE))
  #nodelabels(NodeDepth(tree, T, FALSE))
  expect_equal(c(1,3,2,1,4,4,3,1,2,1,3,1,2,1,3,1,2,1),
               NodeDepth(tree, includeTips = FALSE))
  expect_equal(c(1,2,1,1,3,3,2,1,1,1,2,1,1,1,2,1,1,1),
               NodeDepth(tree, shortest = TRUE, includeTips = FALSE))



  tree <- CollapseNode(BalancedTree(20), c(22:26, 33:35))
  expect_equal(c(4,3,2,1,1,4,1,3,2,1,1), NodeDepth(tree, FALSE, FALSE))
  expect_equal(c(1,2,1,1,1,2,1,2,1,1,1), NodeDepth(tree, TRUE, FALSE))

  tree <- CollapseNode(BalancedTree(20), c(22, 33:35))
  expect_equal(c(4,3,2,1,1,3,2,1,1,4,1,3,2,1,1),
               NodeDepth(tree, FALSE, FALSE))
  expect_equal(c(3,2,1,1,1,2,1,1,1,2,1,2,1,1,1),
               NodeDepth(tree, TRUE, FALSE))

  expect_equal(1L, NodeDepth(CollapseNode(BalancedTree(8), 10:15), FALSE, FALSE))
  expect_equal(1L, NodeDepth(CollapseNode(BalancedTree(8), 10:15), TRUE, FALSE))
})
