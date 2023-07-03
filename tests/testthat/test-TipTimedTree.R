test_that("TipTimedTree() example works", {
  expect_equal(
    TipTimedTree(BalancedTree(6), tipAge = 1:6, minEdge = 2)[["edge.length"]],
    c(5, 2, 3, 2, 3, 2, 2, 3, 2, 3)
  )
  expect_equal(
    TipTimedTree(Postorder(BalancedTree(6)), tipAge = 1:6, minEdge = .5)$edge.l,
    c(0.5, 1.5, 1, 0.5, 0.5, 1.5, 1, 0.5, 0.5, 3.5)
  )
  expect_error(TipTimedTree(BalancedTree(5), 1:10), "one age per leaf")
  expect_error(TipTimedTree(BalancedTree(5), double(0)), "one age per leaf")
  expect_warning(expect_equal(TipTimedTree(SingleTaxonTree(1), 1, 100)$edge.l,
                              100), "does not contain multiple edges")
})
