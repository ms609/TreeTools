test_that("match.multiPhylo()", {
  t1 <- BalancedTree(6)
  t2 <- PectinateTree(7)
  t9 <- StarTree(99)
  trees <- c(t2, t1, StarTree(5), t2, t1)
  
  expect_true(t1 %in% t1)
  expect_false(t9 %in% t1)
  expect_false(t1 %in% t9)
  expect_true(t2 %in% trees)
  expect_false(t9 %in% trees)
  expect_equal(trees %in% t1, c(FALSE, TRUE, FALSE, FALSE, TRUE))
  
  expect_equal(trees[1:2] %in% trees, c(TRUE, TRUE))
  expect_equal(trees %in% trees[1:2], c(TRUE, TRUE, FALSE, TRUE, TRUE))
  
  
  expect_equal(match(t1, t9, nomatch = 123), 123)
  expect_equal(match(t1, t1, nomatch = 123), 1)
  expect_equal(match(t2, c(t1, t2)), 2)
  expect_equal(match(t9, c(t1, t2)), NA_integer_)

  expect_equal(match(c(t1, t2), t2),c(NA, 1))
  expect_equal(match(c(t1, t2), t9), rep(NA_integer_, 2))
  
  expect_equal(match(c(t1, t2), c(t2, t1, t9)), match(1:2, 2:0))
  expect_equal(match(c(t1, t2, t9), c(t2, t1)), match(1:3, 2:1))
  expect_equal(match(c(t1, t2, t1), c(t9, t9)), match(1:3, 5:6))
})
