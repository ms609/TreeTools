nasty <- structure(list(edge = structure(
  c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
    5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
  .Dim = c(12, 2)),
  Nnode = 5L,
  tip.label = letters[1:8]),
  class = "phylo") # Danger: Do not plot!

test_that(".DescendantTips() recurses", {
  nTip <- 19L
  tree <- BalancedTree(nTip)
  edge <- tree$edge
  expect_equal(.DescendantTips(edge[, 1], edge[, 2], nTip, c(23, 32, 35)),
               c(1:3, 11:13, 16:19))
})

test_that("RootOnNode() works", {
  expect_null(RootOnNode(NULL, 1))

  tree <- structure(list(edge = structure(c(6L, 9L, 9L, 7L, 7L, 8L, 8L,
                                            6L, 9L, 2L, 7L, 3L, 8L, 4L, 5L, 1L),
                                          .Dim = c(8L, 2L)),
                         tip.label = c("t3", "t4", "t1", "t2", "t5"),
                         Nnode = 4L), class = "phylo", order = "cladewise")
  
  exp8 <- structure(list(edge = structure(c(6L, 7L, 8L, 8L, 7L, 6L, 9L, 9L, 7L,
                                            8L, 1L, 2L, 3L, 9L, 4L, 5L),
                                          .Dim = c(8L, 2L)),
                         tip.label = c("t3", "t4", "t1", "t2", "t5"),
                         Nnode = 4L), class = "phylo", order = "preorder")
  exp7 <- structure(list(edge = structure(c(6L, 7L, 7L, 6L, 8L, 8L, 9L, 9L, 7L,
                                            1L, 2L, 8L, 3L, 9L, 4L, 5L),
                                          .Dim = c(8L, 2L)),
                         tip.label = c("t3", "t4", "t1", "t2", "t5"),
                         Nnode = 4L), class = "phylo", order = "preorder")
  exp5 <- structure(list(edge = structure(c(6L, 7L, 8L, 9L, 9L, 8L, 7L, 6L, 7L,
                                            8L, 9L, 1L, 2L, 3L, 4L, 5L),
                                          .Dim = c(8L, 2L)),
                         tip.label = c("t3", "t4", "t1", "t2", "t5"),
                         Nnode = 4L), class = "phylo", order = "preorder")
  exp4 <- structure(list(edge = structure(c(6L, 7L, 8L, 9L, 9L, 8L, 7L, 6L, 7L,
                                            8L, 9L, 1L, 2L, 3L, 5L, 4L),
                                          .Dim = c(8L, 2L)),
                         tip.label = c("t3", "t4", "t1", "t2", "t5"),
                         Nnode = 4L), class = "phylo", order = "preorder")
  exp3 <- structure(list(edge = structure(c(6L, 7L, 8L, 8L, 7L, 9L, 9L, 6L, 7L,
                                            8L, 1L, 2L, 9L, 4L, 5L, 3L),
                                          .Dim = c(8L, 2L)),
                         tip.label = c("t3", "t4", "t1", "t2", "t5"),
                         Nnode = 4L), class = "phylo", order = "preorder")
  exp2 <- structure(list(edge = structure(c(6L, 7L, 7L, 8L, 8L, 9L, 9L, 6L, 7L,
                                            1L, 8L, 3L, 9L, 4L, 5L, 2L),
                                          .Dim = c(8L, 2L)),
                         tip.label = c("t3", "t4", "t1", "t2", "t5"),
                         Nnode = 4L), class = "phylo", order = "preorder")
  #t2 <- Preorder(t2)
  expect_true(all.equal(tree, RootOnNode(tree, node = 9L, TRUE)))
  expect_true(all.equal(exp8, RootOnNode(tree, node = 8L, TRUE)))
  expect_true(all.equal(exp7, RootOnNode(tree, node = 7L, TRUE)))
  expect_true(all.equal(tree, RootOnNode(tree, node = 6L, TRUE)))

  expect_true(all.equal(exp5, RootOnNode(tree, node = 5L, TRUE)))
  expect_true(all.equal(exp4, RootOnNode(tree, node = 4L, TRUE)))
  expect_true(all.equal(exp3, RootOnNode(tree, node = 3L, TRUE)))
  expect_true(all.equal(exp2, RootOnNode(tree, node = 2L, TRUE)))
  expect_true(all.equal(tree, RootOnNode(tree, node = 1L, TRUE)))
  expect_true(all.equal(tree, RootOnNode(UnrootTree(tree), node = 1L, TRUE)))

  expb6_8 <- structure(list(edge = structure(c(8, 8, 10, 10, 9, 9, rep(7, 3),
                                               1:2, 4:6, 10, 3, 8:9),
                                             .Dim = c(9L, 2L)),
                            tip.label = paste0("t", 1:6), Nnode = 4L),
                       class = "phylo", order = "preorder")
  expp6_8 <- structure(list(edge = structure(c(7, rep(7:10, each = 2), 1:2, 8,
                                               3, 9, 4, 10, 5:6),
                                             .Dim = c(9L, 2L)),
                            tip.label = paste0("t", 1:6), Nnode = 4L),
                       class = "phylo", order = "preorder")

  expect_true(all.equal(expb6_8, Preorder(RootOnNode(BalancedTree(6L), 8L))))
  expect_true(all.equal(expp6_8, Preorder(RootOnNode(PectinateTree(6L), 8L))))

  expect_true(all.equal(UnrootTree(PectinateTree(6L)),
                        RootOnNode(PectinateTree(6L), 1)))

  TestTip <- function(tr, node, rr) {
    expect_true(all.equal(
      Preorder(ape::root(tr, outgroup = node, resolve.root = rr)),
      RootOnNode(tr, node, rr))
    )
  }
  TestInternal <- function(tr, node, rr) {
    expect_true(all.equal(
      Preorder(ape::root(tr, node = node, resolve.root = rr)),
      RootOnNode(tr, node, rr)))
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

  expect_true(all.equal(PectinateTree(8), RootOnNode(PectinateTree(8), 9L, TRUE)))
  expect_true(all.equal(UnrootTree(PectinateTree(8)),
                        RootOnNode(PectinateTree(8), 9L, FALSE)))
  expect_true(all.equal(urt, RootOnNode(urt, 9L, FALSE)))
  expect_true(all.equal(Preorder(RootTree(urt, letters[1:2])),
                        RootOnNode(urt, 9L, TRUE)))

})

test_that("root_on_node() works", {
  tree <- Preorder(BalancedTree(15))
  edge <- tree$edge
  TipTest <- function(i) {
    tr.rooted <- root_on_node(tree, i)
    expect_equal(SortTree(root(tree, i, resolve.root = TRUE)),
                 SortTree(tr.rooted))
  }
  StaticTest <- function(i) expect_equal(tree, root_on_node(tree, i))
  NodeTest <- function(i) {
    tr.rooted <- root_on_node(tree, i)
    expect_equal(SortTree(root(tree, node = i, resolve.root = TRUE)),
                 SortTree(tr.rooted))
  }
  expect_error(root_on_node(edge, 0))
  for (i in 1:15) TipTest(i)
  StaticTest(16)
  StaticTest(17)
  for (i in 18:23) NodeTest(i)
  StaticTest(24)
  for (i in 24:29) NodeTest(i)
  expect_error(root_on_node(edge, 30))

  tree <- Preorder(root(BalancedTree(15), "t1", resolve.root = TRUE))
  edge <- tree$edge
  expect_error(root_on_node(edge, 0))
  for (i in 1:15) TipTest(i)
  StaticTest(16)
  StaticTest(17)
  for (i in 18:29) NodeTest(i)
  expect_error(root_on_node(edge, 30))
})

test_that("RootOnNode() supports lists of trees", {
  rootOn <- 8L
  expect_equal(structure(list(RootOnNode(as.phylo(1, 5), rootOn),
                              RootOnNode(as.phylo(2, 5), rootOn)),
                         class = "multiPhylo",
                         tip.label = paste0("t", 1:5)),
               RootOnNode(as.phylo(1:2, 5), rootOn))
})

test_that("RootTree() supports star trees", {
  star <- StarTree(8)
  expect_true(all.equal(star, RootTree(star, 5:7)))
  expect_true(all.equal(star$edge, RootTree(star$edge, 5:7)))
})

test_that("RootOnNode() supports nasty node ordering", {
  expect_equal(Preorder(nasty),
               RootOnNode(nasty, 12L, resolveRoot = TRUE))
  expect_equal(RootOnNode(Preorder(nasty), 11L),
               RootOnNode(nasty, 13L))
})

test_that("RootTree() handles null outgroups", {
  bal8 <- BalancedTree(8)
  expect_equal(RootTree(bal8), bal8)
  expect_equal(RootTree(bal8$edge), bal8$edge)
  expect_equal(RootTree(bal8, NULL), bal8)
  expect_equal(RootTree(c(bal8, bal8)), c(bal8, bal8))
  expect_equal(RootTree(list(bal8, bal8)), list(bal8, bal8))
  expect_equal(RootTree(bal8$edge, NULL), bal8$edge)
  expect_equal(RootTree(bal8, character(0)), bal8)
  expect_null(RootTree(NULL))
  expect_null(RootTree(NULL, "tip"))
})

test_that("RootTree() works", {
  bal8 <- BalancedTree(8)
  bal15 <- BalancedTree(15)
  expect_error(RootTree(bal8, 1:8 %in% 0))
  expect_error(RootTree(bal8, "tip_not_there"))
  expect_equal(RootTree(bal8, 5:6), RootTree(bal8, 1:8 %in% 5:6))
  expect_equal(RootTree(bal8$edge, 5:6), RootTree(bal8, 5:6)$edge)
  expect_equal(RootTree(bal15$edge, 9:11), RootTree(bal15, 9:11)$edge)
  expect_equal(RootTree(bal8, 5:6), RootTree(bal8, c("t5", "t6")))
  expect_true(all.equal(
               as.phylo(5518, 8, paste0("t", rev(c(7,8,3,4,1,2,6,5)))),
               RootTree(bal8, "t5")))
  expect_equal(RootTree(bal8, 5), RootTree(bal8, "t5"))
  expect_equal(RootTree(bal8, 5)$edge, RootTree(bal8$edge, 5))
  expect_equal(RootTree(bal8, c("t1", "t2")), RootTree(bal8, c("t4", "t5")))

  expect_equal(structure(list(RootTree(as.phylo(1, 5), "t5"),
                              RootTree(as.phylo(2, 5), "t5")),
                         class = "multiPhylo",
                         tip.label = paste0("t", 1:5)),
               RootTree(as.phylo(1:2, 5), "t5"))

  expect_true(all.equal(read.tree(text = "((b, c), (a, (d, e)));"),
               RootTree(read.tree(text = "(a, ((b, c), (d, e)));"), c(1, 4))))

  pec4 <- PectinateTree(4)
  expect_true(all.equal(pec4, RootTree(pec4, c(1, 3))))
  pec5 <- PectinateTree(5)
  expect_true(all.equal(pec5, RootTree(pec5, c(1, 4))))

  tree <- structure(list(edge = structure(c(7L, 8L, 8L, 7L, 7L, 9L, 9L, 10L,
                                            10L, 8L, 1L, 6L, 2L, 9L, 3L, 10L,
                                            4L, 5L), .Dim = c(9L, 2L)),
                         Nnode = 4L, tip.label = letters[1:6]),
                    class = "phylo", order = "preorder")
  expect_equal(RootTree(tree, 1:5), RootTree(tree, 6))
  # non-contiguous outgroup
  expect_equal(RootTree(tree, 1:4), RootTree(tree, 5:6))
})

test_that("RootTree(fallback)", {
  tangle <- ape::read.tree(text = "((a1, (a2, a3)), ((b4, (b5, b6)), (c7, (c8, c9))));")
  outgroup <- c("a3", "b6", "c7", "c8", "c9")
  fallback <- c("c7", "c8", "b4", "b5", "c9", "a1", "a2")
  toN <- function(x) as.integer(substr(x, 2, 2))

  expect_equal(
    RootTree(tangle, outgroup, fallback),
    RootTree(tangle, c("a3", "b6"))
  )
  expect_equal(
    RootTree(tangle, outgroup, toN(fallback)),
    RootTree(tangle, c("a3", "b6"))
  )
  expect_equal(
    RootTree(tangle, toN(outgroup), fallback),
    RootTree(tangle, c("a3", "b6"))
  )
  expect_equal(
    RootTree(c(tangle, tangle), outgroup, fallback),
    RootTree(c(tangle, tangle), c("a3", "b6"))
  )
})

test_that("RootTree() & UnrootTree() retain edge lengths", {
  bal7 <- BalancedTree(7)
  bal7$edge.length <- 1:12 * 10
  attr(bal7, "order") <- NULL
  expect_equal(RootTree(bal7, 1:4),
               structure(bal7, order = "preorder"))
  expect_equal(RootTree(RootTree(bal7, 1), 1:4),
               structure(bal7, order = "preorder"))
  exp <- structure(bal7, order = "preorder")
  exp$edge.length = 10 * c(1 + 8, 2:7, 0, 9:12)
  expect_equal(RootTree(UnrootTree(bal7), 1:4), exp)
})

test_that("UnrootTree() works", {
  expect_null(UnrootTree(NULL))
  expect_equal(matrix(c(7, 8, 8, 7, 7, 9, 10, 10, 9,
                        8, 1, 2, 3, 9, 10, 4,  5, 6), ncol = 2L),
               UnrootTree(BalancedTree(6))$edge)
  expect_equal(matrix(c(7, 7, 7, 8, 8, 9, 9, 10, 10,
                        1, 2, 8, 3, 9, 4, 10, 5,  6), ncol = 2L),
               UnrootTree(PectinateTree(6))$edge)
  expect_true(all.equal(BalancedTree(2), UnrootTree(BalancedTree(2))))
  expect_true(all.equal(BalancedTree(1), UnrootTree(BalancedTree(1))))
  expect_equal(matrix(c(9, 9, 10, 10, 10, 9, 11, 11, 12, 12, 9,
                        1, 10, 2,  4,  7, 11, 3, 12,  5,  6, 8), ncol = 2L),
               UnrootTree(nasty)$edge)
  # Unrooting when already unrooted
  expect_equal(UnrootTree(PectinateTree(6)),
               UnrootTree(UnrootTree(PectinateTree(6))))

  expList <- list(UnrootTree(as.phylo(1, 5)), UnrootTree(as.phylo(2, 5)))
  expect_equal(expList, UnrootTree(list(as.phylo(1, 5), as.phylo(2, 5))))
  exp <- structure(expList, class = "multiPhylo")
  attr(exp, "tip.label") <- paste0("t", 1:5)
  expect_equal(exp, UnrootTree(as.phylo(1:2, 5)))
})

test_that("CollapseNode() works", {
  tree8  <- read.tree(text="(((a, (b, (c, d))), (e, f)), (g, h));")
  expect_error(CollapseNode(1:5, tree8))
  expect_error(CollapseNode(tree8, 1))
  expect_warning(CollapseNode(tree8, 9L))

  tree <- as.phylo(123, 7)
  tree$edge.length <- 12:1
  expect_true(all.equal(tree, CollapseNode(tree, integer(0))))

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

  expect_true(all.equal(CollapseNode(Preorder(nasty), c(11, 12)),
                        Preorder(CollapseNode(nasty, c(11, 13)))))

  expect_error(CollapseEdge(tree, 9))
})

test_that("CollapseNode() handles node labels", {
  bal6 <- BalancedTree(6)
  startLabels <- paste("Node", 7:11)
  bal6[["node.label"]] <- startLabels
  if (interactive()) {
    plot(bal6, show.node.label = TRUE)
  }
  
  # Collapse a cherry
  expect_equal(
    CollapseNode(bal6, 9)[["node.label"]],
    startLabels[-3]
  )
  
  # Collapse an internal node
  expect_equal(
    CollapseNode(bal6, 8)[["node.label"]],
    startLabels[-2]
  )
  
  # Collapse an internal node
  expect_equal(
    CollapseEdge(bal6, 1)[["node.label"]],
    startLabels[-2]
  )
  
  # case = 3 -> y is bound on an internal edge
  expect_equal(
    CollapseNode(bal6, c(8, 11))[["node.label"]],
    startLabels[-c(2, 5)]
  )
  
  expect_equal(AddTipEverywhere(bal6)[[1]][["node.label"]],
               AddTip(bal6, where = 1)[["node.label"]])
})

test_that("MakeTreeBinary() edge lengths", {
  tree <- ape::read.tree(text = "(a:1, b:2, c:3, d:4);")
  # Not yet supported; ensure they are removed
  expect_null(MakeTreeBinary(tree)$edge.length)
})

test_that("Binarification is uniform", {
  set.seed(0)
  Test <- function(tree, nTree, nSamples = 200L, ape = FALSE) {
    counts <- table(replicate64(nSamples, as.TreeNumber(MakeTreeBinary(tree))))
    expect_equal(nTree, length(counts))
    expect_gt(chisq.test(counts)$p.value, 0.001)
  }

  Test(CollapseNode(PectinateTree(6), 8:9), NUnrooted(4))
  Test(CollapseNode(PectinateTree(6), 9:10), NRooted(4))
  Test(CollapseNode(PectinateTree(6), c(8, 10)), NUnrooted(3) * NRooted(3))
  Test(CollapseNode(BalancedTree(8), c(10:12)), NUnrooted(5))
  Test(CollapseNode(BalancedTree(7), c(10, 13)), NRooted(3) * NRooted(3))

  bal7 <- BalancedTree(7)
  bal7[["node.label"]] <- paste("Node", 8:13)
  expect_true(all.equal(bal7, MakeTreeBinary(bal7)))
  expect_true(all.equal(list(bal7, bal7), MakeTreeBinary(list(bal7, bal7))))
  expect_true(all.equal(
    structure(list(bal7, bal7), class = "multiPhylo"),
    MakeTreeBinary(structure(list(bal7, bal7), class = "multiPhylo"))))

  set.seed(1)
  binNodes <- MakeTreeBinary(CollapseNode(bal7, 9:10))[["node.label"]]
  expect_equal(binNodes[!is.na(binNodes)], bal7[["node.label"]][-(2:3)])
  
})

test_that("LeafLabelInterchange() fails", {
  skip_if(Sys.getenv("USING_ASAN") != "")
  expect_error(LeafLabelInterchange(BalancedTree(4), 5)) # n too many
})

test_that("LeafLabelInterchange() works", {
  expect_equal(PectinateTree(40), LeafLabelInterchange(PectinateTree(40), 1))
  expect_true(all.equal(BalancedTree(2),
                        LeafLabelInterchange(BalancedTree(2), 2)))
  expect_equal(rev(BalancedTree(2)$tip),
                   LeafLabelInterchange(BalancedTree(2), 2)$tip)

  abcd <- letters[1:4]
  sapply(1 + seq_len(100), function(i) {
    # Check all perturbations
    set.seed(i)
    expect_false(any(abcd == LeafLabelInterchange(BalancedTree(abcd), 4)$tip))

    # Check lots of sizes
    expect_equal(i, sum(paste0("t", 1:128) !=
                        LeafLabelInterchange(BalancedTree(128), i)$tip))
  })
})
