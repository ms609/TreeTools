context("splits.cpp")

test_that("Bad data results in error", {
  expect_error(cpp_edge_to_splits(matrix(1, 3, 3), 3))
  expect_error(cpp_edge_to_splits(matrix(1, 4e6, 2), 3))
  expect_error(cpp_edge_to_splits(matrix(1, 10, 2), 0))
  expect_error(cpp_edge_to_splits(matrix(1, 10, 2), -10))
  expect_error(cpp_edge_to_splits(matrix(1, 2, 2), 6)) # Impossible anyway...
})
