test_that("Treeness() fails safely", {
  expect_error(Treeness(list(1, 2)), "list of .phylo. objects")
  expect_warning(expect_null(Treeness(BalancedTree(8))), "edge lengths")
  expect_null(Treeness(NULL))
})

test_that("Treeness() succeeds", {
  internal8 <- c(1, 2, 5, 8, 9, 12)
  in8 <- 1:14 %in% internal8
  tr0 <- BalancedTree(8, !in8)
  tr1 <- BalancedTree(8, in8)
  expect_equal(Treeness(tr1), 1)
  expect_equal(Stemminess(tr0), 0)
  expect_equal(Treeness(BalancedTree(8, ifelse(in8, 1, 11))),
               1 * sum(in8) / sum(ifelse(in8, 1, 11)))
  expect_equal(Stemminess(c(tr0, tr1)), 0:1)
  expect_equal(Stemminess(list(tr0, tr1)), 0:1)
  expect_error(Stemminess(list(tr0, tr1, "tr1")), "phylo")
})
