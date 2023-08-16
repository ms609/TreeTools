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
