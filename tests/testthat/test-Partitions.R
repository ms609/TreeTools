context("Partitions.R")

test_that("Split combination", {
  tree1 <- ape::read.tree(text='((a, b), (c, (d, e)));')
  splits1 <- as.Splits(tree1)
  tree2 <- ape::read.tree(text='(a, (b, (c, (d, e))));')
  splits12 <- c(splits1, as.Splits(tree2))
  tree3 <- ape::read.tree(text='(a, (b, (c, (e, d))));')
  tree4 <- ape::read.tree(text='(a, (b, (c, (d, f))));')
  tree5 <- ape::read.tree(text='(a, (b, (c, (d, (e, f)))));')

  expect_equal(4L, length(splits12))
  expect_equal(c(FALSE, FALSE, TRUE, TRUE), as.logical(duplicated(splits12)))
  expect_error(c(splits1, as.Splits(tree3)))
  expect_error(c(splits1, as.Splits(tree4)))
  expect_error(c(splits1, as.Splits(tree5)))
  expect_equal(c(28L, 24L, 28L, 24L),
               as.integer(c(splits1,
                            as.Splits(RenumberTips(tree3, letters[1:5])))[, 1]))
  expect_equal(2L, length(unique(c(splits1, as.Splits(tree2)))))
  expect_equal(2L, length(unique(c(splits1,
                                   as.Splits(c(3, 7), tipLabels=letters[1:5])))))

  # TODO: Fully test splits with large (>32 tip) trees
})
