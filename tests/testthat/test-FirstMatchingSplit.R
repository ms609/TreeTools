test_that("FirstMatchingSplit() fails gracefully", {
  bal13 <- BalancedTree(13)
  pec13 <- PectinateTree(13)
  expect_error(FirstMatchingSplit(bal13, pec13),
               "Splits")
})

test_that("FirstMatchingSplit() works", {
  bal13 <- as.Splits(BalancedTree(13))
  pec13 <- as.Splits(PectinateTree(13))
  expect_equal(
    FirstMatchingSplit(bal13, pec13),
    which(bal13 %in% pec13)[[1]]
  )
  
  expect_equal(FirstMatchingSplit(bal13, as.Splits(StarTree(13))), 0)
  expect_equal(FirstMatchingSplit(as.Splits(StarTree(13)), bal13, nomatch = NA),
               NA_integer_)
})
