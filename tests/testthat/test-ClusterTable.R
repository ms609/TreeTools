test_that("ClusterTable fails gracefully", {
  bigTree <- PectinateTree(2^14 + 1)
  expect_error(
    as.ClusterTable(bigTree),
    "Tree has too many leaves. Contact the .TreeTools. maintainer."
  )
})

test_that("ClusterTable class behaves", {
  tree <- RootTree(BalancedTree(6), 1)
  ct <- as.ClusterTable(tree)
  expect_equal(matrix(c(0, 2, rep(1, 4), 0, 3, 3:6), 6),
               as.matrix.ClusterTable(ct))
  expect_equal(ClusterTable_decode(ct), 6:1)
  
  expect_equal(capture.output(ct),
               "ClusterTable on 6 leaves: t1 .. t6")
  
  expect_equal(
    capture.output(summary(ct)),
    c("ClusterTable on 6 leaves:",
      " 123456",
      " .**...",
      " ***...",
      " ****..",
      " *****.",
      " ******",
      " 1: t6  2: t5  3: t4  4: t3  5: t2  6: t1 ")
  )
  
  ClusterSummary <- function(...) capture.output(summary(as.ClusterTable(...)))
  t1..8 <- paste0("t", 1:8)
  t8..1 <- rev(t1..8)
  expect_equal(
    ClusterSummary(BalancedTree(t8..1),
                   tipLabels = t1..8),
    ClusterSummary(BalancedTree(t1..8))
  )
  
  byList <- as.ClusterTable(list(BalancedTree(t8..1), PectinateTree(t8..1)),
                            tipLabels = t1..8)
  listBy <- list(as.ClusterTable(BalancedTree(t1..8)),
                 as.ClusterTable(PectinateTree(t1..8)))
  expect_equal(capture.output(summary(byList[[2]])),
               capture.output(summary(listBy[[2]])))
  
})

test_that("ClusterTable with multiple trees", {
  tree1 <- ape::read.tree(text = "(A, (B, (C, (D, E))));");
  tree2 <- ape::read.tree(text = "(E, (B, (D, (C, A))));");
  ct <- as.ClusterTable(c(tree1, tree2))
  expect_equal(
    capture.output(print(ct)),
    capture.output(print(list(as.ClusterTable(tree1),
                              as.ClusterTable(tree2, TipLabels(tree1)))))
  )
  expect_equal(
    capture.output(summary(ct)),
    capture.output(summary(list(as.ClusterTable(tree1),
                                as.ClusterTable(tree2, TipLabels(tree1)))))
  )
  
  trees <- ape::read.tree(text = 
                            "(A, (B, (C, (D, E))));(E, (B, (D, (C, A))));");
  expect_equal(
    capture.output(print(as.ClusterTable(trees))),
    capture.output(print(list(as.ClusterTable(trees[[1]]),
                              as.ClusterTable(trees[[2]], TipLabels(tree1)))))
  )
})

test_that("Attributes are correct", {
  t6 <- as.ClusterTable(BalancedTree(6))
  t7 <- as.ClusterTable(PectinateTree(7))
  t8 <- as.ClusterTable(BalancedTree(8))
  s8 <- StarTree(8)
  expect_equal(3, NSplits(t6))
  expect_equal(4:5, NSplits(list(t7, t8)))
  
  expect_equal(6, NTip(t6))
  expect_equal(7:8, NTip(list(t7, t8)))
  
  #TODO test TipLabels, SplitsInBalancedTree
})

test_that("ClusterTable with multiple trees", {
  tree1 <- ape::read.tree(text = "(A, (B, (C, (D, E))));");
  tree2 <- ape::read.tree(text = "(E, (B, (D, (C, A))));");
  ct <- as.ClusterTable(c(tree1, tree2))
  expect_equal(
    capture.output(print(ct)),
    capture.output(print(list(as.ClusterTable(tree1),
                              as.ClusterTable(tree2, TipLabels(tree1)))))
  )
  expect_equal(
    capture.output(summary(ct)),
    capture.output(summary(list(as.ClusterTable(tree1),
                                as.ClusterTable(tree2, TipLabels(tree1)))))
  )
  
  trees <- ape::read.tree(text = 
                            "(A, (B, (C, (D, E))));(E, (B, (D, (C, A))));");
  expect_equal(
    capture.output(print(as.ClusterTable(trees))),
    capture.output(print(list(as.ClusterTable(trees[[1]]),
                              as.ClusterTable(trees[[2]], TipLabels(tree1)))))
  )
})

test_that("ClusterTable with complex trees", {
  skip_if_not_installed("TreeDist", "2.9.2.9000")
  library("TreeDist")
  
  # Test exposes failures in C++ - constexpr not playing nicely with Rcpp
  # Specifically if replacing
  # const int16 L_COL = int16(0);
  # const int16 R_COL = int16(1);
  # const int16 X_COLS = int16(2);
  # hard-coding, using enum all fails.
  tr1 <- structure(list(
    edge = structure(c(8L, 8L, 9L, 10L, 10L, 9L, 11L, 11L, 8L, 12L, 12L,
                       1L, 9L, 10L, 2L, 3L, 11L, 4L, 5L, 12L, 6L, 7L),
                     dim = c(11L, 2L)),
    Nnode = 5L, tip.label = c("t1", "t2", "t3", "t4", "t5", "t6", "t7")),
    class = "phylo", order = "preorder")
  tr2 <- structure(list(
    edge = structure(c(8L, 9L, 10L, 10L, 9L, 11L, 11L, 8L, 12L, 12L, 8L,
                       9L, 10L, 1L, 2L, 11L, 3L, 4L, 12L, 5L, 6L, 7L),
                     dim = c(11L, 2L)),
    Nnode = 5L, tip.label = c("t1", "t2", "t3", "t4", "t5", "t6", "t7")),
    class = "phylo", order = "preorder")
  t4 <- list(a = tr1, b = tr2, c = tr1, d = tr2)
  r4 <- RootTree(t4, 1)
  
  expect_equal(as.numeric(RobinsonFoulds(r4)), c(8, 0, 8, 8, 0, 8))
})
