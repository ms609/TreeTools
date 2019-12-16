context("TreeNumber.R")

test_that("as.phylo.numeric", {
  expect_equal(as.phylo(0:2, 6, letters[1:6])[[1]],
               PectinateTree(letters[c(1, 3:6, 2)]))
})
