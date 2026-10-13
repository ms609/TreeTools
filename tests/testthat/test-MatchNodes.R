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

test_that("MatchNodes() handles differing tip sets", {
  tree <- BalancedTree(8)
  # Clades must contain exactly the same tips
  sub <- DropTip(tree, "t8")
  expect_equal(MatchNodes(sub, tree), c(NA, 10:12, NA, 14))
  expect_equal(MatchNodes(tree, sub), c(NA, 9:11, NA, 13, NA))
  # ... but MatchEdges() disregards tips of `table` absent from `x`
  expect_equal(MatchEdges(sub, tree), 1:12)
  expect_equal(MatchEdges(tree, sub, nomatch = 0L),
               c(1:7, 0, 9:11, 0, 12, 0))
  
  # Tip order is immaterial
  shuffled <- RenumberTips(tree, rev(TipLabels(tree)))
  expect_equal(MatchNodes(tree, shuffled, tips = TRUE), c(8:1, 9:15))
})

test_that("MatchNodes() handles nodes with a single child", {
  # Unary nodes above clade (t1, t2), above tip t3 and at the root
  tree <- read.tree(text = "((((t1, t2)), (t3)), t4);")
  ref <- read.tree(text = "(((t1, t2), t3), t4);")
  expect_equal(MatchNodes(ref, tree), c(5, 6, 7))
  expect_equal(MatchNodes(tree, ref), c(5, 6, 7, 7, 3))
  expect_equal(MatchNodes(tree, ref, tips = TRUE),
               c(1:4, 5, 6, 7, 7, 3))
  # First matching edge is reported
  expect_equal(MatchEdges(ref, tree), c(1, 2, 4, 5, 6, 8))
  rootUnary <- read.tree(text = "((((t1, t2), t3), t4));")
  expect_equal(MatchNodes(rootUnary, ref), c(5, 5, 6, 7))
  expect_equal(MatchNodes(ref, rootUnary), c(6, 7, 8))
})

test_that("MatchNodes() matches reference implementation", {
  set.seed(161)
  for (i in 1:40) {
    x <- RandomTree(sample(4:30, 1), root = runif(1) < 0.5)
    table <- switch(sample(4, 1),
                    RandomTree(TipLabels(x), root = TRUE),
                    RootTree(x, sample(TipLabels(x), 1)),
                    CollapseNode(x, NTip(x) + 2),
                    DropTip(x, sample(TipLabels(x), 2)))
    table <- RenumberTips(table, sample(TipLabels(table)))
    expect_equal(
      MatchNodes(x, table, tips = TRUE),
      .MatchNodesR(x[["edge"]], table[["edge"]], TipLabels(x),
                   TipLabels(table), tips = TRUE)
    )
    expect_equal(
      MatchEdges(x, table),
      .MatchEdgesR(x[["edge"]], table[["edge"]], TipLabels(x),
                   TipLabels(table))
    )
  }
})

test_that("MatchNodes() falls back for duplicated labels", {
  dup <- BalancedTree(c("a", "b", "c", "c"))
  bal4 <- BalancedTree(c("a", "b", "c", "d"))
  expect_equal(MatchNodes(dup, bal4),
               .MatchNodesR(dup[["edge"]], bal4[["edge"]],
                            TipLabels(dup), TipLabels(bal4)))
  expect_equal(MatchEdges(dup, bal4),
               .MatchEdgesR(dup[["edge"]], bal4[["edge"]],
                            TipLabels(dup), TipLabels(bal4)))
})

test_that("MatchNodes() rejects malformed trees", {
  bal4 <- BalancedTree(4)
  twoParents <- bal4
  twoParents[["edge"]][2, 2] <- twoParents[["edge"]][1, 2]
  expect_error(MatchNodes(bal4, twoParents), "two parents")
  expect_error(MatchNodes(twoParents, bal4), "two parents")
  expect_error(match_clades(matrix(1:3, 1), bal4[["edge"]], 1:4, FALSE),
               "two columns")
  expect_error(match_clades(bal4[["edge"]] - 1L, bal4[["edge"]], 1:4, FALSE),
               "positive")
})
