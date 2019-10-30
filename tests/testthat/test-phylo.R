context('phylo.R')

test_that('EverywhereTips correct', {
  backbone <- PectinateTree(5)

  expect_equal(7L, length(AddTipEverywhere(backbone, includeRoot = FALSE)))
  expect_equal(9L, length(AddTipEverywhere(backbone, includeRoot = TRUE)))

})
