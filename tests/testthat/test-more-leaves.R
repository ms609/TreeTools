test_that("big trees can be compared", {
  nTip <- 10000
  bal <- unroot(BalancedTree(nTip))
  balSplits <- as.Splits(bal)
  expect_equal(dim(balSplits), c(NTip(bal) - 3, ceiling(NTip(bal) / 8)))
  expect_equal(Preorder(splits_to_edge(balSplits, nTip)), Preorder(bal$edge))
})

test_that("Big tree distances can be computed", {
  skip_if_not_installed("TreeDist")
  bal <- BalancedTree(10000)
  pec <- PectinateTree(10000)
  TreeDist::TreeDistance(bal, pec)
})
