nasty <- structure(list(edge = structure(
  c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
    5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
  .Dim = c(12, 2)),
  Nnode = 5L,
  tip.label = letters[1:8]),
  class = 'phylo') # Danger: Do not plot!


context("Tree rearrangements")

test_that("RootOnNode works", {

  tree <- structure(list(edge = structure(c(6L, 9L, 9L, 7L, 7L, 8L, 8L,
                                            6L, 9L, 2L, 7L, 3L, 8L, 4L, 5L, 1L),
                                          .Dim = c(8L, 2L)),
                         tip.label = c("t3", "t4", "t1", "t2", "t5"), Nnode = 4L),
                    class = "phylo", order = "cladewise")

  exp8 <- structure(list(edge = structure(c(6L, 7L, 8L, 8L, 7L, 6L, 9L, 9L, 7L, 8L, 1L, 2L, 3L, 9L, 4L, 5L), .Dim = c(8L, 2L)), tip.label = c("t3", "t4", "t1", "t2", "t5"), Nnode = 4L), class = "phylo", order = "preorder")
  exp7 <- structure(list(edge = structure(c(6L, 7L, 7L, 6L, 8L, 8L, 9L, 9L, 7L, 1L, 2L, 8L, 3L, 9L, 4L, 5L), .Dim = c(8L, 2L)), tip.label = c("t3", "t4", "t1", "t2", "t5"), Nnode = 4L), class = "phylo", order = "preorder")
  exp5 <- structure(list(edge = structure(c(6L, 7L, 8L, 9L, 9L, 8L, 7L, 6L, 7L, 8L, 9L, 1L, 2L, 3L, 4L, 5L), .Dim = c(8L, 2L)), tip.label = c("t3", "t4", "t1", "t2", "t5"), Nnode = 4L), class = "phylo", order = "preorder")
  exp4 <- structure(list(edge = structure(c(6L, 7L, 8L, 9L, 9L, 8L, 7L, 6L, 7L, 8L, 9L, 1L, 2L, 3L, 5L, 4L), .Dim = c(8L, 2L)), tip.label = c("t3", "t4", "t1", "t2", "t5"), Nnode = 4L), class = "phylo", order = "preorder")
  exp3 <- structure(list(edge = structure(c(6L, 7L, 8L, 8L, 7L, 9L, 9L, 6L, 7L, 8L, 1L, 2L, 9L, 4L, 5L, 3L), .Dim = c(8L, 2L)), tip.label = c("t3", "t4", "t1", "t2", "t5"), Nnode = 4L), class = "phylo", order = "preorder")
  exp2 <- structure(list(edge = structure(c(6L, 7L, 7L, 8L, 8L, 9L, 9L, 6L, 7L, 1L, 8L, 3L, 9L, 4L, 5L, 2L), .Dim = c(8L, 2L)), tip.label = c("t3", "t4", "t1", "t2", "t5"), Nnode = 4L), class = "phylo", order = "preorder")
  #t2 <- Preorder(t2)
  expect_equal(tree, RootOnNode(tree, node = 9L, TRUE))
  expect_equal(exp8, RootOnNode(tree, node = 8L, TRUE))
  expect_equal(exp7, RootOnNode(tree, node = 7L, TRUE))
  expect_equal(tree, RootOnNode(tree, node = 6L, TRUE))

  expect_equal(exp5, RootOnNode(tree, node = 5L, TRUE))
  expect_equal(exp4, RootOnNode(tree, node = 4L, TRUE))
  expect_equal(exp3, RootOnNode(tree, node = 3L, TRUE))
  expect_equal(exp2, RootOnNode(tree, node = 2L, TRUE))
  expect_equal(tree, RootOnNode(tree, node = 1L, TRUE))
  expect_equal(tree, RootOnNode(unroot(tree), node = 1L, TRUE))


  TestTip <- function (tr, node, rr) {
    expect_equal(Preorder(ape::root(tr, outgroup = node, resolve.root = rr)),
                 RootOnNode(tr, node, rr))
  }
  TestInternal <- function (tr, node, rr) {
    expect_equal(Preorder(ape::root(tr, node = node, resolve.root = rr)),
                 RootOnNode(tr, node, rr))
  }

  TestInternal(PectinateTree(8), 14, TRUE)
  TestInternal(PectinateTree(8), 14, FALSE)
  TestTip(PectinateTree(8), 4L, TRUE)
  TestTip(PectinateTree(8), 4L, FALSE)

  TestInternal(BalancedTree(8L), 11L, TRUE)
  TestInternal(BalancedTree(8L), 11L, FALSE)

  urt <- UnrootedTreeWithShape(3, 8, letters[1:8])
  TestInternal(urt, 12, TRUE)
  TestInternal(urt, 12, FALSE)
  TestTip(urt, 4L, TRUE)
  TestTip(urt, 4L, FALSE)

  # Children of an unresolved root
  TestInternal(urt, 10, TRUE)
  TestTip(urt, 1, TRUE)

  expect_equal(PectinateTree(8), RootOnNode(PectinateTree(8), 9L, TRUE))
  expect_equal(unroot(PectinateTree(8)), RootOnNode(PectinateTree(8), 9L, FALSE))
  expect_equal(urt, RootOnNode(urt, 9L, FALSE))
  expect_equal(Preorder(EnforceOutgroup(urt, letters[1:2])),
               RootOnNode(urt, 9L, TRUE))
})

test_that("RootOnNode supports nasty node ordering", {
  expect_equal(Preorder(nasty),
               RootOnNode(nasty, 12L, resolveRoot = TRUE))
  expect_equal(RootOnNode(Preorder(nasty), 11L),
               RootOnNode(nasty, 13L))
})

test_that("CollapseNodes works", {
  tree8  <- read.tree(text="(((a, (b, (c, d))), (e, f)), (g, h));")
  expect_error(CollapseNode(1:5, tree8))
  expect_error(CollapseNode(tree8, 1))
  expect_warning(CollapseNode(tree8, 9L))

  tree <- as.phylo(123, 7)
  tree$edge.length <- 12:1
  expect_equal(tree, CollapseNode(tree, integer(0)))

  exp1213 <- matrix(c(8, 8, 9, 10, 10, 9, 11, 11, 11, 11,
                      1, 9, 10, 2, 4, 11, 3, 7, 6, 5), ncol=2)
  no1213 <- CollapseNode(tree, c(12, 13))
  expect_equal(exp1213, no1213$edge)

  el <- tree$edge.length
  expect_equal(no1213$edge.length, c(el[1:6],
                                     el[7] + c(c(el[8] + el[9:10]), el[11]),
                                     el[12]))

  no11 <- CollapseEdge(tree, c(7, 8))
  expect_equal(exp1213, no11$edge)

  expect_equal(CollapseNode(Preorder(nasty), c(11, 12)),
               Preorder(CollapseNode(nasty, c(11, 13))))
})

test_that("DropTip works", {
  bal8 <- BalancedTree(8)
  expect_null(DropTip(bal8, 1:8))
  expect_warning(expect_equal(bal8, DropTip(bal8, -1)))
  expect_warning(expect_equal(bal8, DropTip(bal8, 99)))
  expect_warning(expect_equal(bal8, DropTip(bal8, 'MissingTip')))
  expect_error(DropTip(bal8, list('Invalid format')))

  expect_equal(DropTip(bal8, 7:8), DropTip(bal8, 15L))
  expect_equal(ape::drop.tip(bal8, 6:8), DropTip(bal8, 6:8))
  expect_equal(ape::drop.tip(bal8, c(3, 5, 7)), DropTip(bal8, c(3, 5, 7)))

  expect_equal(DropTip(Preorder(nasty), c(1, 3)),
               Preorder(DropTip(nasty, c(1, 3))))


  bigTree <- RandomTree(1284)
  set.seed(1284)
  bigTip <- sample(1:1284, 608)
  expect_equal(ape::drop.tip(bigTree, bigTip), DropTip(bigTree, bigTip))
  #microbenchmark(ape::drop.tip(bigTree, bigTip), DropTip(bigTree, bigTip), times = 25)
  #profvis(replicate(25, DropTip(bigTree, bigTip)), interval = 0.005)
})
