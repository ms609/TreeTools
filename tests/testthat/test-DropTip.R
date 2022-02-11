test_that("keep_tip() works", {
  expect_error(keep_tip(matrix(1:4, 4, 1), FALSE))
})
