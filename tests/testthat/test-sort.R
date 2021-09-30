test_that("Trees are sorted", {
  trees <- as.phylo(0:8, 8)
  expect_true(trees[[1]] == trees[[1]])
  expect_false(trees[[1]] == trees[[2]])
  expect_true(trees[[1]] < trees[[2]] || trees[[1]] > trees[[2]])
})
