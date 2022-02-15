test_that("Bad data results in error", {
  skip_if(Sys.getenv("USING_ASAN") != "")
  expect_error(cpp_edge_to_splits(matrix(1, 3, 3), 1:3, 3))
  expect_error(cpp_edge_to_splits(matrix(1, 10, 2), 1:10, 0))
  expect_error(cpp_edge_to_splits(matrix(1, 10, 2), 1:10, -10))
  expect_error(cpp_edge_to_splits(matrix(1, 10, 2), 1:9, 5))
  expect_error(cpp_edge_to_splits(matrix(1, 2, 2), 1:2, 6)) # Impossible anyway...
})
