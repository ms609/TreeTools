test_that("KeptVerts() works", {
  bal12 <- BalancedTree(12)
  Y <- TRUE
  N <- FALSE
  expect_equal(KeptVerts(bal12, 1:12 %in% 1:9),
               c(1:12 %in% 1:9, rep(Y, 6), N, Y, Y, N, N))
  expect_equal(KeptVerts(bal12, 1:12 %in% 7:12),
               c(1:12 %in% 7:12, rep(N, 6), rep(Y, 5)))
  expect_equal(KeptVerts(bal12, 1:12 %in% 4:6),
               c(1:12 %in% 4:6, rep(N, 4), rep(Y, 2), rep(N, 5)))
  expect_equal(KeptVerts(bal12, 1:12 %in% c(1, 2, 6, 12)),
               c(1:12 %in% c(1, 2, 6, 12), Y, Y, N, Y, N, N, N, N, N, N, N))
})
