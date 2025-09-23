test_that("descendant_edges() handles errors", {
  expect_error(descendant_edges(1:2, 1:3, 1:2),
               "`parent` and `child` must be the same length")
  expect_error(descendant_tips(1:2, 1:3, 1:2),
               "`parent` and `child` must be the same length")
  expect_error(descendant_edges(1:2, 1:2, 1:3),
               "`postorder` must list each edge once")
  expect_error(descendant_tips(1:2, 1:2, 1:3),
               "`postorder` must list each edge once")
})

test_that("DescendantEdges() works", {
  pec5 <- UnrootTree(PectinateTree(5))
  V <- TRUE
  x <- FALSE
  answer <- matrix(c(V, x, x, x, x, x, x,
                     x, V, x, x, x, x, x,
                     x, x, V, V, V, V, V,
                     x, x, x, V, x, x, x,
                     x, x, x, x, V, V, V,
                     x, x, x, x, x, V, x,
                     x, x, x, x, x, x, V), 7, 7, byrow = T)
  expect_equal(
    DescendantEdges(edge = NULL, pec5$edge[, 1], pec5$edge[, 2]),
    answer)
  expect_equal(
    apply(DescendantEdges(node = 0, pec5$edge[, 1], pec5$edge[, 2]), 1, which),
    list(1:7, 4:7, 6:7)
  )
  expect_equal(
    apply(DescendantEdges(node = 8:7, pec5$edge[, 1], pec5$edge[, 2]), 1, which),
    list(6:7, 4:7)
  )
})

test_that("DescendantTips() works", {
  tree <- as.phylo(0, 6)
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  expect_equal(DescendantTips(parent, child, 5), 1:6 %in% c(2, 6))
  expect_equal(
    DescendantTips(parent, child),
    t(vapply(1:10, function(x) DescendantTips(parent, child, edge = x),
             logical(6)))
  )
  
  polytomies <- CollapseNode(BalancedTree(9), c(12, 13, 16))
  if (interactive()) {
    plot(polytomies)
    edgelabels()
    nodelabels()
  }
  pol <- polytomies[["edge"]]
  expect_equal(which(DescendantTips(pol[, 1], pol[, 2], edge = 5)), 4:5)
  expect_equal(which(DescendantTips(pol[, 1], pol[, 2], edge = 1)), 1:5)
  expect_equal(which(DescendantTips(pol[, 1], pol[, 2], edge = 8)), 6:9)
  expect_equal(which(DescendantTips(pol[, 1], pol[, 2], edge = 10)), 7)
  
  expect_equal(which(DescendantTips(pol[, 1], pol[, 2], node = 5)), 5)
  expect_equal(which(DescendantTips(pol[, 1], pol[, 2], node = 11)), 1:5)
  expect_equal(which(DescendantTips(pol[, 1], pol[, 2], node = 0)[11, ]), 1:5)
})

test_that("DescendantTips() handles postorder", {
  post6 <- Postorder(BalancedTree(6))
  if (interactive()) {
    oPar <- par(mar = rep(0.5, 4), cex = 0.9)
    on.exit(par(oPar))
    plot(post6)
    nodelabels()
    edgelabels()
    tiplabels()
  }
  parent <- post6[["edge"]][, 1]
  child <- post6[["edge"]][, 2]
  expect_equal(DescendantTips(parent, child, edge = 7), 1:6 %in% 1:2)
  expect_equal(DescendantTips(parent, child, edge = 9), 1:6 %in% 4:6)
  # edge = NULL is handled by .AllDescendantEdges
  expect_equal(DescendantTips(parent, child),
               t(vapply(list(5, 4, 4:5, 6, 2, 1, 1:2, 3, 4:6, 1:3),
                        function (table) 1:6 %in% table,
                        logical(6)))
  )
})
