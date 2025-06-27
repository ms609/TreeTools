test_that("Reweight()", {
  mat <- rbind(a = c(c1 = 0, c2 = 2, c3 = 0), b = c(0, 2, 0), c = c(1, 3, 0),
               d = c(1, 3, 0))
  dat <- MatrixToPhyDat(mat)
  
  expect_error(Reweight(mat, letters[1:3]),
               "must be a numeric vector")
  expect_error(Reweight(mat, 1:2),
               "Length of `weights` \\(2\\) must match .* characters \\(3\\)")
  
  expect_equal(
    Reweight(mat, c(1, 2, 0)),
    `colnames<-`(mat[, c(1, 2, 2)], c("c1_1", "c2_1", "c2_2"))
  )
  # Equivalently:
  Reweight(dat, c("3" = 0, "2" = 2))
})