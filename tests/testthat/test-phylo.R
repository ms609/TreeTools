context('phylo.R')

nasty <- structure(list(edge = structure(
  c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
    5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
  .Dim = c(12, 2)),
  Nnode = 5L,
  tip.label = letters[1:8]),
  class = 'phylo') # Danger: Do not plot!

test_that('AddTipEverywhere correct', {
  backbone <- PectinateTree(5)

  expect_equal(7L, length(AddTipEverywhere(backbone, includeRoot = FALSE)))
  expect_equal(7L, length(AddTipEverywhere(unroot(backbone), includeRoot = FALSE)))
  expect_equal(9L, length(AddTipEverywhere(backbone, includeRoot = TRUE)))
  expect_equal(5L, length(AddTipEverywhere(CollapseNode(backbone, 7:9))))

})

test_that('AddTipEverywhere handles nasty tree', {
  expect_equal(AddTipEverywhere(Preorder(nasty)),
               lapply(AddTipEverywhere(nasty), Preorder))
})

test_that('ListAncestors works', {
  edge <- nasty$edge
  expect_equal(c(10L, 12L), ListAncestors(edge[, 1], edge[, 2], 11))
})

test_that("CladeSizes works", {
  #plot(Preorder(nasty)); nodelabels(c(12, 10, 13, 11, 9)); tiplabels(1:8)
  #edgelabels(c(2, 3, 6, 4, 8, 10, 9, 5, 7, 1, 12, 11))

  expect_equal(c(3, 8 + 4, 3 + 1, 7 + 3, 2),
               CladeSizes(nasty, internal = TRUE, 13:9))
  expect_equal(c(3, 8, 3, 7, 2), CladeSizes(nasty, internal = FALSE, 13:9))
})
