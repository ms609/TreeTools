test_that("ReadMrBayes() fails gracefully", {
  expect_error(MrBayesTrees("notA-file", burninFrac = 1),
               "Cannot find notA-file\\.run1")
})
