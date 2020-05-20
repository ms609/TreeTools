context("int_to_tree.cpp")

test_that("Failures are graceful", {
  expect_error(num_to_parent(10, 1))
  expect_error(num_to_parent(10, -1))
  expect_error(edge_to_num(1:10, 1:11, 6))
  expect_error(edge_to_num(1:10, 1:10, 5))
})
