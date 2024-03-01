test_that("Inf cons runs", {
  trees <- as.phylo(0:5, 8)
  if (interactive()) {
    plot(trees[[1]])
    nodelabels(c(0, 0:5))
  }
  Consensus(trees)
})
