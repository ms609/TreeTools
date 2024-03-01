test_that("Inf cons runs", {
  trees <- as.phylo(0:5, 8)
  if (interactive()) {
    plot(trees[[1]])
    nodelabels(c(0, 0:5))
  }
  counts <- count_splits(trees)
  expect_equal(sum(counts$count), NSplits(trees[[1]]) * length(trees))
  expect_equal(ncol(counts$splits), NTip(trees[[1]]))
  splits <- as.Splits(counts$splits)
  expect_false(any(duplicated(splits)))
  
  skip_if(TRUE)
  cons <- Consensus(trees)
  expect_equal(mode(cons), "list")
  expect_equal(length(cons), 4)
  dput(cons)
})
