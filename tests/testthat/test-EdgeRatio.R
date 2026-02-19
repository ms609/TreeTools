test_that("EdgeRatio() works", {
  expect_warning(expect_equal(EdgeRatio(BalancedTree(5)), NA_real_),
                 "[Ee]dge lengths not spec")
  
  ex2 <- BalancedTree(4)
  ex2[["edge.length"]] <- c(1, 2, 2, 1, 2, 2)
  expect_equal(EdgeRatio(ex2),
               structure(8/2, external = 8, internal = 2, total = 10))
})
