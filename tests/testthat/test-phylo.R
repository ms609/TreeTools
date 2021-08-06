nasty <- structure(list(edge = structure(
  c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
    5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
  .Dim = c(12, 2)),
  Nnode = 5L,
  tip.label = letters[1:8]),
  class = 'phylo') # Danger: Do not plot!

test_that('AddTipEverywhere() correct', {
  backbone <- PectinateTree(5)

  expect_equal(7L, length(AddTipEverywhere(backbone, includeRoot = FALSE)))
  expect_equal(7L, length(AddTipEverywhere(UnrootTree(backbone),
                                           includeRoot = FALSE)))
  expect_equal(9L, length(AddTipEverywhere(backbone, includeRoot = TRUE)))
  expect_equal(5L, length(AddTipEverywhere(CollapseNode(backbone, 7:9))))

})

test_that("AddTip() at root", {
  expect_equal(PectinateTree(8),
               AddTip(DropTip(PectinateTree(8), 1), where = 0, label = 't1'))
  AddTip(nasty, 1L)
  AddTip(nasty, 12L)
})

test_that("AddTip() with tip name", {
  bal8 <- BalancedTree(8)
  expect_error(AddTip(bal8, 'invalid tip'))
  expect_equal(AddTip(bal8, 1L), AddTip(bal8, 't1'))
  expect_equal(AddTip(nasty, 1L), AddTip(nasty, 'a'))
})

test_that("AddTip() with edge lengths", {
  pec8 <- PectinateTree(8)
  pec8$edge.length <- rep(1L, 14L)

  # case = 1 -> y is bound on the root of x
  expect_equal(c(0.6, rep(1, 14), 2),
               AddTip(pec8, 0, edgeLength = 2, lengthBelow = 0.6)$edge.length)
  # case = 2 -> y is bound on a tip of x
  expect_equal(c(rep(1, 2), 0.3, 0.7, 2, rep(1, 11)),
               AddTip(pec8, 2, edgeLength = 2, lengthBelow = 0.7)$edge.length)
  # case = 3 -> y is bound on a node of x
  expect_equal(c(rep(1, 11), 0.3, 2, 0.7, 1, 1),
               AddTip(pec8, 15, edgeLength = 2, lengthBelow = 0.7)$edge.length)
  expect_equal(c(rep(1, 11), 0.5, 0, 0.5, 1, 1),
               AddTip(pec8, 15)$edge.length)
  expect_equal(c(rep(1, 11), 0, 2, 1, 1, 1),
               AddTip(pec8, 15, edgeLength = 2, lengthBelow = 1)$edge.length)
  expect_equal(c(rep(1, 11), -0.6, 2, 1.6, 1, 1),
               AddTip(pec8, 15, edgeLength = 2, lengthBelow = 1.6)$edge.length)
})

test_that('AddTipEverywhere() handles nasty tree', {
  added <- AddTipEverywhere(nasty)
  lapply(added, function (tr) expect_true(all(tr$edge > 0)))
  expect_equal(AddTipEverywhere(Preorder(nasty)),
               lapply(added, Preorder))
})

test_that('AddTipEverywhere() with tiny trees', {
  added <- AddTipEverywhere(StarTree(2))
  lapply(added, function (tr) expect_true(all(tr$edge > 0)))
  expect_equal(2, length(added))
  expect_equal(3, length(AddTipEverywhere(StarTree(2), include = TRUE)))

  expect_equal(list(PectinateTree(c('t1', 'New tip'))),
               AddTipEverywhere(StarTree(1)))
  expect_equal(list(SingleTaxonTree('New tip')),
               AddTipEverywhere(structure(list(tip.label = character(0)),
                                          class = 'phylo')))
})

test_that("Subtree() works", {
  expect_error(Subtree(BalancedTree(8), 10)) # Nodes must be in preorder
  t4 <- Subtree(Preorder(BalancedTree(8)), 10)
  expect_equal(BalancedTree(4), t4)
  expect_equal(BalancedTree(4), Subtree(t4, 5))
  expect_equal(SingleTaxonTree('t1'), Subtree(t4, 1))
})


test_that('ListAncestors() works', {
  edge <- nasty$edge
  expect_equal(c(10L, 12L), ListAncestors(edge[, 1], edge[, 2], 11))
  tr <- PectinateTree(4)
  expect_equal(list(5, 6:5, 7:5, 7:5, integer(0), 5, 6:5),
               ListAncestors(tr$edge[, 1], tr$edge[, 2]))
  expect_equal(integer(0), ListAncestors(tr$edge[, 1], tr$edge[, 2], 5))
})

test_that("CladeSizes() works", {
  #plot(Preorder(nasty)); nodelabels(c(12, 10, 13, 11, 9)); tiplabels(1:8)
  #edgelabels(c(2, 3, 6, 4, 8, 10, 9, 5, 7, 1, 12, 11))

  expect_equal(c(3, 8 + 4, 3 + 1, 7 + 3, 2),
               CladeSizes(nasty, internal = TRUE, 13:9))
  expect_equal(c(3, 8, 3, 7, 2), CladeSizes(nasty, internal = FALSE, 13:9))

  # Misspecification:
  expect_warning(expect_equal(CladeSizes(BalancedTree(7), internal = FALSE),
                              CladeSizes(BalancedTree(7), internal = 8:9)))
})
