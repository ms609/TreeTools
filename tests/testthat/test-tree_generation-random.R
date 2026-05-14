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
  ape_counts <- table(vapply64(lapply(rep(5, nSamples), ape::rtree),
                               as.TreeNumber, 1))
  expect_lt(chisq.test(ape_counts)$p.value, 0.005)

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

test_that("RandomTree(root = TRUE) samples uniformly from rooted binary trees", {
  # There are (2n-3)!! labelled rooted binary trees on n tips:
  #   n=4: 15 trees; n=5: 105 trees.
  # Prior to fix, root_on_node sampled nodes (not edges), giving only 3 distinct
  # topologies at n=4.  Regression test: chi-sq GoF against uniform.
  nReps <- 1e4
  canon <- function(tr) paste(ape::write.tree(Preorder(tr)))

  set.seed(1)
  trees4 <- lapply(rep(4L, nReps), RandomTree, root = TRUE)
  counts4 <- table(vapply(trees4, canon, character(1L)))
  expect_gt(length(counts4), 10L)   # all 15 topologies reachable
  expect_gt(chisq.test(counts4)$p.value, 0.001)

  set.seed(2)
  trees5 <- lapply(rep(5L, nReps), RandomTree, root = TRUE)
  counts5 <- table(vapply(trees5, canon, character(1L)))
  expect_gt(length(counts5), 50L)   # well above 0 of 105 topologies reachable
  expect_gt(chisq.test(counts5)$p.value, 0.001)
})
