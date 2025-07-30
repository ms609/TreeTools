test_that("LongBranch() fails safely", {
  expect_error(LongBranch(list(1, 2)), "list of .phylo. objects")
  expect_warning(expect_null(LongBranch(BalancedTree(8))), "edge lengths")
  expect_null(LongBranch(NULL))
})

test_that("LongBranch() succeeds", {
  evenLength <- BalancedTree(8, rep(1, 14))
  expect_equal(unname(LongBranch(evenLength)), rep(0, 8))
  tree <- BalancedTree(8, lengths = c(rep(2, 4), 5:7, rep(2, 4), rep(1, 3)))
  lb <- LongBranch(tree)
  pat <- ape::cophenetic.phylo(tree)
  expect_equal(lb, (colSums(pat) / (NTip(tree) - 1) /
                      mean(pat[upper.tri(pat)])) - 1)
  expect_equal(LongBranch(c(tree, tree)), list(lb, lb))
  expect_equal(LongBranch(c(evenLength, tree)),
               LongBranch(list(evenLength, tree)))
})
