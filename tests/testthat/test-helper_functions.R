context("helper_functions.R")

test_that("UnshiftTree() works", {
  t1 <- as.phylo(1, 8)
  t2..9 <- as.phylo(2:9, 8)
  t1..9 <- as.phylo(1:9, 8)
  expect_equivalent(t1..9, UnshiftTree(t1, t2..9))
  expect_equivalent(t1..9, UnshiftTree(t1, lapply(t2..9, I)))
  expect_equivalent(as.phylo(1:2, 8), UnshiftTree(t1, as.phylo(2, 8)))
})
