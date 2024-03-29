test_that("MatchNodes() works", {
  bal8 <- BalancedTree(8)
  bal16 <- BalancedTree(16)
  if (interactive()) {
    oPar <- par(mfrow = c(2, 1), mar = rep(1, 4), xpd = NA)
    on.exit(par(oPar))
    plot(bal8)
    edgelabels()
    nodelabels()
    par(cex = 0.75)
    plot(bal16)
    edgelabels()
    nodelabels()
  }
  expect_equal(MatchEdges(bal8, bal16), 2:15)
  e16 <- bal16[["edge"]]
  expect_equal(.Edge(bal8$edge), bal8[["edge"]])
  l16 <- structure(list(e16[, 1], e16[, 2]), tip.label = TipLabels(bal16))
  expect_equal(.Edge(l16), e16)
  
  expect_equal(MatchEdges(l16, bal8), c(NA, 1:14, rep(NA, 15)))
  expect_equal(MatchEdges(bal16$edge, bal8$edge), c(NA, 1:14, rep(NA, 15)))
  expect_equal(MatchNodes(bal8, bal16, tips = FALSE), 18:24)
  expect_equal(MatchNodes(bal8, bal16, tips = TRUE), c(1:8, 18:24))
  expect_equal(MatchNodes(BalancedTree(16), bal8, tips = TRUE, nomatch = -1),
               c(1:8, rep(-1, 8), -1, 9:15, rep(-1, 7)))
  
  bal8rn <- bal8
  bal8rn[["tip.label"]][7:8] <- letters[7:8]
  expect_equal(MatchNodes(bal8rn, bal8, tips = TRUE),
               c(1:6, rep(NA, 2), NA, 10:12, NA, 14, NA))
  
  table <- RootTree(bal8, 1)
  if (interactive()) {
    par(cex = 0.9)
    plot(bal8)
    edgelabels()
    nodelabels()
    plot(table)
    edgelabels()
    nodelabels()
  }
  expect_equal(MatchNodes(bal8, table),
               c(9, rep(NA_integer_, 2), 12:15))
  expect_equal(MatchNodes(bal8, table, nomatch = -1),
               c(9, rep(-1, 2), 12:15))
  
  
  if (interactive()) {
    plot(bal8)
    edgelabels()
    nodelabels()
    plot(Postorder(bal8))
    edgelabels()
    nodelabels()
  }
  expect_equal(MatchNodes(bal8, Postorder(bal8)),
               c(9:15))
})

test_that(".UpdateNodeLabel() works", {
  .Node <- function(n) paste("Node", n)
  bal4 <- BalancedTree(4)
  bal4[["node.label"]] <- .Node(5:7)
  bal8 <- BalancedTree(8)
  bal8[["node.label"]] <- .Node(9:15)
  expect_equal(.UpdateNodeLabel(bal4, bal8), .Node(10:12))
  expect_equal(.UpdateNodeLabel(bal8, bal4), c(NA, .Node(5:7), rep(NA, 3)))
  
  expect_equal(.UpdateNodeLabel(bal4$edge, bal8, newTips = TipLabels(4)),
               .Node(10:12))
  expect_equal(.UpdateNodeLabel(bal8$edge, bal4, newTips = TipLabels(8)),
               c(NA, .Node(5:7), rep(NA, 3)))
  
})