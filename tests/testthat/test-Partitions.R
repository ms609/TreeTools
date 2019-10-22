context("Partitions.R")

test_that("Split combination", {
  tree1 <- ape::read.tree(text='((a, b), (c, (d, e)));')
  splits1 <- as.Splits(tree1)
  tree2 <- ape::read.tree(text='(a, (b, (c, (d, e))));')
  tree3 <- ape::read.tree(text='(a, (b, (c, (e, d))));')
  tree4 <- ape::read.tree(text='(a, (b, (c, (d, f))));')
  tree5 <- ape::read.tree(text='(a, (b, (c, (d, (e, f)))));')

  expect_equal(4L, length(c(splits1, as.Splits(tree2))))
  expect_error(c(splits1, as.Splits(tree3)))
  expect_error(c(splits1, as.Splits(tree4)))
  expect_error(c(splits1, as.Splits(tree5)))
  expect_equal(c(28L, 24L, 28L, 24L),
               as.integer(c(splits1,
                            as.Splits(RenumberTips(tree3, letters[1:5])))[, 1]))
  expect_equal(2L, length(unique(c(splits1, as.Splits(tree2)))))
  expect_equal(2L, length(unique(c(splits1,
                                   as.Splits(c(3, 7), tipLabels=letters[1:5])))))

})

test_that("UniqueSplits works", {
  suppressWarnings(RNGversion("3.5.0")) # Until we can require R3.6.0
  set.seed(1)
  splits6 <- Tree2Splits(ape::rtree(6, br=NULL))
  expect_equal(c('8'=FALSE, '10'=FALSE, '11'=FALSE), UniqueSplits(splits6)['t4', ])
  expect_equal(!splits6, UniqueSplits(cbind(!splits6, splits6), TRUE))

})

test_that("Large splits don't cause memory issues", {
  splits5000 <- cbind(c(rep(TRUE, 2), rep(FALSE, 4998)),
                      c(rep(FALSE, 2), rep(TRUE, 4998)),
                      c(rep(TRUE, 2500), rep(FALSE, 2500)),
                      c(rep(TRUE, 2500), rep(FALSE, 2500)),
                      c(rep(TRUE, 2500), rep(FALSE, 2500)),
                      c(rep(FALSE, 2500), rep(TRUE, 2500))
  )
  expect_equal(c(5000, 2), dim(UniqueSplits(splits5000)))
})
