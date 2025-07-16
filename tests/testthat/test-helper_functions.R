test_that("UnshiftTree() works", {
  t1 <- as.phylo(1, 8)
  t2..9 <- setNames(as.phylo(2:3, 8), letters[2:3])
  t1..9 <- setNames(as.phylo(1:3, 8), c("", letters[2:3]))
  attr(t1..9, "tip.label") <- NULL
  expect_true(all.equal(UnshiftTree(t1, t2..9), t1..9))
  expect_equal(unclass(t1..9), UnshiftTree(t1, unclass(t2..9)))
  expectation <- as.phylo(1:2, 8)
  attr(expectation, "tip.label") <- NULL
  expect_equal(expectation, UnshiftTree(t1, as.phylo(2, 8)))
})
