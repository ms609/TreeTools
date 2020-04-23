context("tree_properties.R")

test_that('EdgeAncestry works', {
  tree <- BalancedTree(10)
  edge <- tree$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  expect_equal(c(10, 14, 16), which(EdgeAncestry(18, parent, child)))
  expect_equal(1, which(EdgeAncestry(2, parent, child)))
  expect_equal(integer(0), which(EdgeAncestry(1, parent, child)))
  expect_equal(logical(18), EdgeAncestry(10, parent, child))
})

test_that('Root node can be found', {
  rooted <- BalancedTree(8)
  postorder <- Postorder(rooted)
  unrooted <- unroot(rooted)
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

test_that('NodeOrder works', {
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

  tree8 <- ape::rtree(8L)
  expect_true(TreeIsRooted(tree8))
  expect_false(TreeIsRooted(ape::unroot(tree8)))

  expect_equal(5L, NPartitions(8))
  expect_equal(5L, NPartitions(tree8))
  expect_equal(5L, NPartitions(ape::unroot(tree8)))

  expect_equal(0L, NPartitions(ape::read.tree(text='(A, ((B, C)));')))

  expect_equal(c(5L, 5L), NPartitions(list(tree8, ape::unroot(tree8))))
  expect_equal(c(5L, 5L), NPartitions(c(8, 8)))
  expect_error(NPartitions('not a tree'))
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
  #plot(PectinateTree(20))
  #nodelabels(NodeDepth(PectinateTree(20), F))
  expect_equal(c(rep(0, 20), 1:10, 9:1), NodeDepth(PectinateTree(20)))

  #plot(BalancedTree(20))
  #nodelabels(NodeDepth(BalancedTree(20), F))
  expect_equal(c(5,4,3,1,2,1,3,1,2,1,4,3,1,2,1,3,1,2,1),
               NodeDepth(BalancedTree(20L), FALSE))

  #tree <- CollapseNode(BalancedTree(20), c(22:26, 33:35))
  #plot(tree)
  #nodelabels(NodeDepth(tree, FALSE))
  expect_equal(c(4,3,1,2,1,4,1,3,1,2,1), NodeDepth(tree, FALSE))

  #tree <- CollapseNode(BalancedTree(20), c(22, 33:35))
  #plot(tree)
  #nodelabels(NodeDepth(tree, FALSE))
  expect_equal(c(4,3,1,2,1,3,1,2,1,4,1,3,1,2,1),
               NodeDepth(tree, FALSE))

  expect_equal(1L, NodeDepth(CollapseNode(BalancedTree(8), 10:15), FALSE))
})
