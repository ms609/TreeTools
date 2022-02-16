test_that("Preorder survives", {
  edg <- cbind(rep(9, 8), 1:8)
  expect_equal(Postorder(edg), edg)
})
