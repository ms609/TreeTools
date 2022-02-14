test_that("KeptPaths() works", {
  tree <- BalancedTree(7)
  tree$edge.length <- 1:12
  paths <- PathLengths(tree)
  keptTips <- c(2, 4, 6)
  kept <- KeptVerts(tree, 1:7 %in% keptTips)
  expect_equal(paths[KeptPaths(paths, kept), "length"],
               c(1 + 2 + 4,
                 2 + 4,
                 1 + 5 + 7,
                 5 + 7,
                 8 + 9 + 11,
                 1))
  expected <- PathLengths(KeepTip(tree, keptTips))
  expected[, 1:2] <- which(kept)[unlist(expected[, 1:2])]
  
  KeptPaths(PathLengths(tree, TRUE), kept)
})
  