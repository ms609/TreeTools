test_that("SortTree() works", {
  expect_equal(SortTree(as.phylo(10, 6))$edge,
               matrix(c(7:10, 10, 9, 11, 11, 8, 7:10, 3:4, 11, 1:2, 5:6), 10))
  expect_warning(
    expect_equal(
      SortTree(as.phylo(10, 6), how = "unrecognized"),
      SortTree(as.phylo(10, 6))
    )
  )
  
  start <- makeNodeLabel(as.phylo(10, 6))
  if(interactive()) {
    oPar <- par(mfrow = c(2, 1), mar = rep(0.5, 4), cex = 0.9)
    on.exit(par(oPar))
    plot(start, show.node.label = TRUE, xpd = NA)
    nodelabels()
    plot(SortTree(start), show.node.label = TRUE, xpd = NA)
    nodelabels()
  }
  expect_equal(SortTree(start)[["node.label"]], paste0("Node", c(1:3, 5, 4)))
  
  tipSorted <- SortTree(as.phylo(10, 6), how = "tip")
  children <- tipSorted[["edge"]][, 2]
  expect_equal(children[children < 7],
               order(TipLabels(tipSorted))[c(1, 2, 6, 3, 4, 5)])
  
  tipSorted <- SortTree(as.phylo(10, 6), how = "tip",
                        order = paste0("t", c(6, 4, 7, 2)))
  children <- tipSorted[["edge"]][, 2]
  expect_equal(children[children < 7],
               order(TipLabels(tipSorted))[c(6, 2, 4, 3, 5, 1)])
  
  
  unrooted5 <- UnrootTree(PectinateTree(5))
  unrooted5$edge.length <- 1:7
  sorted5 <- SortTree(unrooted5)
  expect_equal(sorted5$edge,
               matrix(c(6, 7, 8, 8, 7, 6, 6,
                        7, 8, order(sorted5[["tip.label"]], decreasing = TRUE)),
                      ncol = 2))
  expect_equal(sorted5[["edge.length"]], c(3, 5, 7, 6, 4, 2, 1))
  expect_equal("cladewise", attr(SortTree(PectinateTree(5)), "order"))
  expect_equal(SortTree(sorted5, "tip",
                        paste0("t", c(1, 2, 4, 5, 3)))[["edge.length"]],
               c(1:3, 5:7, 4))

  
  skip_if_not_installed("vdiffr", "1.0")
  vdiffr::expect_doppelganger("sorted-tree", function() {
    par(mar = rep(0, 4))
    plot(SortTree(UnrootTree(CollapseNode(BalancedTree(17), c(28, 31, 32)))))
    edgelabels(adj = c(1, 1/2))
  })
})


test_that("SortTree.multiPhylo()", {
  t1 <- as.phylo(123, 12)
  t2 <- as.phylo(921, 12)

  expect_identical(structure(list(SortTree(t1), SortTree(t2)),
                             class = "multiPhylo"),
                   SortTree(c(t1, t2)))
})
