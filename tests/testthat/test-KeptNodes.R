test_that("KeptVerts() works", {
  bal12 <- BalancedTree(12)
  expect_error(KeptVerts(bal12, 1:9), "must be a logical vector")
  Y <- TRUE
  N <- FALSE
  expect_error(KeptVerts(matrix(FALSE, 4, 2), c(Y, Y, N, Y)),
               "no applicable method")
  expect_error(KeptVerts(matrix(0, 4, 3), c(Y, Y, N, Y)),
               "edge matrix of a `phylo` obj")
  expect_error(KeptVerts(1:8, c(Y, Y, N, Y)),
               "edge matrix of a `phylo` obj")
  
  expect_equal(KeptVerts(bal12, 1:12 %in% 1:9),
               c(1:12 %in% 1:9, rep(Y, 6), N, Y, Y, N, N))
  expect_equal(KeptVerts(bal12, 1:12 %in% 7:12),
               c(1:12 %in% 7:12, rep(N, 6), rep(Y, 5)))
  expect_equal(KeptVerts(bal12, 1:12 %in% 4:6),
               c(1:12 %in% 4:6, rep(N, 4), rep(Y, 2), rep(N, 5)))
  expect_equal(KeptVerts(bal12, 1:12 %in% c(1, 2, 6, 12)),
               c(1:12 %in% c(1, 2, 6, 12), Y, Y, N, Y, N, N, N, N, N, N, N))
})
