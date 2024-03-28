nasty <- structure(list(edge = structure(
  c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
    5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
  .Dim = c(12, 2)),
  Nnode = 5L,
  tip.label = letters[1:8]),
  class = "phylo") # Danger: Do not plot!

test_that("Subtree() works", {
  expect_error(Subtree(read.tree(text = "((a,b),(c,d));", 7)),
               " be in preorder")
  t4 <- Subtree(Preorder(BalancedTree(8)), 10)
  expect_true(all.equal(BalancedTree(4), t4))
  expect_true(all.equal(BalancedTree(4), Subtree(t4, 5)))
  expect_true(all.equal(SingleTaxonTree("t1"), Subtree(t4, 1)))
})

test_that("Subtree() handles node labels", {
  bal8 <- ape::makeNodeLabel(RootTree(BalancedTree(8), 1), prefix = "Node ")
  expect_equal(Subtree(bal8, 5 + 8)[["node.label"]],
               paste("Node", 5:7))
})


test_that("ListAncestors() works", {
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
