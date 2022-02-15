nasty <- structure(list(edge = structure(
  c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
    5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
  .Dim = c(12, 2)),
  Nnode = 5L,
  tip.label = letters[1:8]),
  class = 'phylo') # Danger: Do not plot!

test_that("AllDescendantEdges() works", {
  pec5 <- UnrootTree(PectinateTree(5))
  V <- TRUE; x <- FALSE
  expect_equal(matrix(c(V, x, x, x, x, x, x,
                        x, V, x, x, x, x, x,
                        x, x, V, V, V, V, V,
                        x, x, x, V, x, x, x,
                        x, x, x, x, V, V, V,
                        x, x, x, x, x, V, x,
                        x, x, x, x, x, x, V), 7, 7, byrow = T),
               DescendantEdges(edge = NULL, pec5$edge[, 1], pec5$edge[, 2]))
})

test_that('EdgeAncestry() works', {
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
  expect_equal(ignore_attr = TRUE, c(2L, rep(3L, 6)),
               as.integer(NodeOrder(BalancedTree(8), internalOnly = TRUE)))
  expect_equal(ignore_attr = TRUE, c(rep(1L, 8), 2L, rep(3L, 6)),
               as.integer(NodeOrder(BalancedTree(8), internalOnly = FALSE)))
  expect_equal(ignore_attr = TRUE, rep(2L, 7),
               as.integer(NodeOrder(BalancedTree(8), FALSE)))
  tree <- CollapseNode(BalancedTree(8), 12:15)
  expect_equal(5:3, as.integer(NodeOrder(tree$edge, internalOnly = TRUE)))
  expect_equal(list(NodeOrder(tree), NodeOrder(tree)), NodeOrder(c(tree, tree)))
  expect_equal(NodeOrder(c(tree, tree)), NodeOrder(list(tree, tree)))

  expect_equal(c(3, 2, 2), NDescendants(UnrootTree(PectinateTree(5))))
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
  expect_equal(2L, NSplits('((a, b), (c, (d, e)));'))
  expect_equal(0L, NSplits('a'))
  expect_null(NSplits(NULL))
  expect_equal(0L, NSplits(letters[1:2]))
  expect_equal(3L, NSplits(letters[1:6]))
  expect_error(NPartitions(raw(1)))
})

test_that("NTip() works", {
  Test <- function(n) {
    tr <- BalancedTree(n)
    pec <- PectinateTree(n)
    expect_identical(NTip(tr), n)
    expect_identical(NTip(tr$edge), n)
    expect_identical(NTip(Postorder(tr$edge)), n)
    expect_identical(NTip(list(tr)), n)
    expect_identical(rep(n, 2L), NTip(list(tr, tr)))
    expect_identical(rep(n, 3L),
                     NTip(structure(list(tr, tr, pec), class = 'multiPhylo')))
    expect_null(NTip(matrix('', n, 2)))
  }
  Test(1L)
  Test(8L)
  Test(64L)
  Test(2000L)
  expect_identical(48L, NTip(Lobo.phy))
})

test_that("NSplits() works", {
  expect_equal(NSplits(5L), NSplits(LETTERS[1:5]))
  expect_equal(3L, NSplits(structure(matrix(0, 6, 2), class = 'ClusterTable')))
})

test_that("MRCA() works", {
  bal7 <- BalancedTree(7)
  allAnc <- AllAncestors(bal7$edge[, 1], bal7$edge[, 2])
  expect_equal(9, MRCA(1, 4, allAnc))
  expect_equal(8, MRCA(1, 6, allAnc))
  expect_equal(8, MRCA(1, 7, allAnc))
  expect_equal(1, MRCA(1, 1, allAnc))
  expect_equal(9, MRCA(1, 11, allAnc))
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
  expect_equal(list(c(4,3,2,1,1,1,2,1,1,1,3,2,1,1,1,2,1,1,1)),
               NodeDepth(list(tree), shortest = TRUE, includeTips = FALSE))

  tree <- UnrootTree(tree)
  expect_equal(c(4,3,2,1,1,3,2,1,1,4,3,2,1,1,3,2,1,1),
               NodeDepth(tree, includeTips = FALSE))
  expect_equal(c(3,2,1,1,1,2,1,1,1,3,2,1,1,1,2,1,1,1),
               NodeDepth(tree, shortest = TRUE, includeTips = FALSE))

  tree <- UnrootTree(RootTree(BalancedTree(20), 't10'))
  #plot(tree); nodelabels(adj = 2, bg = 'yellow')
  #nodelabels(NodeDepth(tree, F, FALSE))
  #nodelabels(NodeDepth(tree, T, FALSE))
  expect_equal(c(1,3,4,3,2,1,1,4,3,2,1,1,3,2,1,1,2,1),
               NodeDepth(tree, includeTips = FALSE))
  expect_equal(c(1,2,3,2,1,1,1,3,2,1,1,1,2,1,1,1,1,1),
               NodeDepth(tree, shortest = TRUE, includeTips = FALSE))



  tree <- CollapseNode(BalancedTree(20), c(22:26, 33:35))
  expect_equal(c(4,3,2,1,1,4,1,3,2,1,1), NodeDepth(tree, FALSE, FALSE))
  expect_equal(c(1,2,1,1,1,2,1,2,1,1,1), NodeDepth(tree, TRUE, FALSE))

  tree <- CollapseNode(BalancedTree(20), c(22, 33:35))
  expect_equal(c(4,3,2,1,1,3,2,1,1,4,1,3,2,1,1),
               NodeDepth(tree, FALSE, FALSE))
  expect_equal(c(3,2,1,1,1,2,1,1,1,2,1,2,1,1,1),
               NodeDepth(tree, TRUE, FALSE))

  expect_equal(1L, NodeDepth(StarTree(8), FALSE, FALSE))
  expect_equal(1L, NodeDepth(StarTree(8), TRUE, FALSE))
})

test_that('SplitsInBinaryTree() works', {
  expect_identical(5L, SplitsInBinaryTree(8))
  expect_identical(5L, SplitsInBinaryTree(as.phylo(0, 8)))
  expect_identical(rep(5L, 4), SplitsInBinaryTree(as.phylo(0:3, 8)))
  expect_null(SplitsInBinaryTree(NULL))
})
