test_that("FirstMatchingSplit() fails gracefully", {
  bal13 <- BalancedTree(13)
  pec13 <- PectinateTree(13)
  expect_error(FirstMatchingSplit(bal13, raw(13)),
               "Splits")
})

test_that("FirstMatchingSplit() works", {
  bal13 <- as.Splits(BalancedTree(13))
  pec13 <- as.Splits(PectinateTree(13))
  
  expect_equal(
    FirstMatchingSplit(BalancedTree(13), PectinateTree(13), return = "both"),
    FirstMatchingSplit(bal13, pec13, return = "both")
  )
  
  
  firstMatch <- which(bal13 %in% pec13)[[1]]
  expect_equal(FirstMatchingSplit(bal13, pec13), firstMatch)
  
  expect_equal(
    FirstMatchingSplit(bal13, pec13, return = "both"),
    c(firstMatch, match(bal13[[firstMatch]], pec13))
  )
  expect_equal(FirstMatchingSplit(bal13, as.Splits(StarTree(13))), 0)
  expect_equal(FirstMatchingSplit(as.Splits(StarTree(13)), bal13, nomatch = NA),
               NA_integer_)
  
  # Check robustness to label order
  ren13 <- RenumberTips(pec13, TipLabels(13:1))[[-1]]
  # Split 1 is t1, t2, t3; not t11, t12, t13
  # So first match is split 2 (t1:4).
  expect_equal(FirstMatchingSplit(ren13, bal13, return = "both"), c(2, 2))
})
