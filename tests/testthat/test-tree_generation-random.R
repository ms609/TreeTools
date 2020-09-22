context('tree_generation.R: RandomTree()')

test_that("Random trees are set by R seed", {

  set.seed(0)
  tr1 <- RandomTree(123)
  tr2 <- RandomTree(123)
  set.seed(0)
  tr3 <- RandomTree(123)
  tr4 <- RandomTree(123)
  expect_identical(tr1, tr3)
  expect_identical(tr2, tr4)
  expect_false(identical(tr1, tr2))

  set.seed(1)
  tr5 <- RandomTree(123)
  set.seed(1)
  tr6 <- RandomTree(123)
  expect_identical(tr5, tr6)
  expect_false(identical(tr1, tr5))

})

test_that("Random trees drawn from uniform distribution", {
  # NB We're also testing sapply64 and vapply64 in this chunk
  expect_error(.RandomParent(1))
  expect_error(.RandomParent(0))
  expect_error(.RandomParent(-1))

  nSamples <- 100
  # Ape's trees are not uniformly distributed:
  counts <- table(vapply64(lapply(rep(5, nSamples), ape::rtree),
                           as.TreeNumber, 1))
  expect_lt(chisq.test(counts)$p.value, 0.001)

  # Our trees are:
  counts <- table(sapply64(lapply(rep(5, nSamples), RandomTree), as.TreeNumber))
  expect_gt(chisq.test(counts)$p.value, 0.001)

  if (FALSE) { # Takes many seconds - but worth checking manually?
    counts <- table(sapply64(lapply(rep(8, 200000), RandomTree), as.TreeNumber))
    expect_gt(chisq.test(counts)$p.value, 0.001)
  }

  expect_false(is.rooted(RandomTree(6, root = FALSE)))
  expect_true(is.rooted(RandomTree(6, root = TRUE)))
})
