test_that("SortTree() works", {
  expect_equal(SortTree(as.phylo(10, 6))$edge,
               matrix(c(7:10, 10, 9, 11, 11, 8, 7:10, 3:4, 11, 1:2, 5:6), 10))
  unrooted5 <- UnrootTree(PectinateTree(5))
  unrooted5$edge.length <- 1:7
  expect_warning(sorted5 <- SortTree(unrooted5)) #49
  expect_equal(sorted5$edge,
               matrix(c(6, 7, 8, 8, 7, 6, 6,
                        7, 8, order(sorted5[["tip.label"]], decreasing = TRUE)),
                      ncol = 2))
  expect_null(sorted5[["edge.length"]])
  expect_equal('cladewise', attr(SortTree(PectinateTree(5)), 'order'))
  
  skip_if_not_installed("vdiffr", "1.0")
  vdiffr::expect_doppelganger('sorted-tree', function() {
    par(mar = rep(0, 4))
    plot(SortTree(UnrootTree(CollapseNode(BalancedTree(17), c(28, 31, 32)))))
    edgelabels(adj = c(1, 1/2))
  })
})

test_that("SortTree.multiPhylo()", {
  t1 <- as.phylo(123, 12)
  t2 <- as.phylo(921, 12)

  expect_identical(structure(list(SortTree(t1), SortTree(t2)),
                             class = 'multiPhylo'),
                   SortTree(c(t1, t2)))
})
