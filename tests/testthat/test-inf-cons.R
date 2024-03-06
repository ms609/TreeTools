test_that("Inf cons runs", {
  trees <- as.phylo(0:5, 8)

  counts <- count_splits(trees)
  all_splits <- as.Splits(counts$splits)
  dup <- duplicated(all_splits)
  
  expect_equal(sum(counts$count[!dup]), NSplits(trees[[1]]) * length(trees))
  expect_equal(ncol(counts$splits), NTip(trees[[1]]))
  
  splits <- all_splits[[!dup]]
  nTrees <- length(trees)
  maj_splits <- splits[[counts[["count"]][!dup] > nTrees / 2]]
  cons_splits <- as.Splits(Consensus(trees, p = 0.5, inf = FALSE))
  
  expect_equal(length(maj_splits), length(cons_splits))
  expect_true(all(maj_splits %in% cons_splits))
})
