test_that("MatchNodes() works", {
  bal8 <- BalancedTree(8)
  bal16 <- BalancedTree(16)
  table <- RootTree(bal8, 1)
  if (interactive()) {
    oPar <- par(mfrow = c(2, 1), mar = rep(1, 4))
    plot(bal8)
    edgelabels()
    nodelabels()
    par(cex = 0.75)
    plot(bal16)
    edgelabels()
    nodelabels()
    par(oPar)
  }
  expect_equal(MatchEdges(bal8, bal16), 2:15)
  expect_equal(MatchEdges(bal16, bal8),
               c(NA, 1:14, rep(NA, 15)))
  expect_equal(MatchNodes(bal8, bal16), c(1:8, 18:24))
  expect_equal(MatchNodes(BalancedTree(16), bal8), 
               c(NA, 9:15, rep(NA, 7)))
  
  
  expect_equal(MatchNodes(bal8, table, nomatch = -1),
               c(-1, -1, -1, 12:15))
  expect_equal(MatchNodes(bal8, table, nomatch = -1),
               c(-1, -1, -1, 12:15))
  expect_equal(MatchNodes(bal8, Postorder(bal8)),
               c(-1, -1, -1, 12:15))
})