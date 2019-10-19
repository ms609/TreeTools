context("tree_properties.R")

test_that('Rooting and partition counting', {
  set.seed(0)
  
  expect_false(TreeIsRooted(ape::read.tree(text='(a, b, c);')))
  expect_true(TreeIsRooted(ape::read.tree(text='(a, (b, c));')))
  
  tree8 <- ape::rtree(8L)
  expect_true(TreeIsRooted(tree8))
  expect_false(TreeIsRooted(ape::unroot(tree8)))
  
  expect_equal(5L, NPartitions(8))
  expect_equal(5L, NPartitions(tree8))
  expect_equal(5L, NPartitions(ape::unroot(tree8)))
  
  expect_equal(0L, NPartitions(ape::read.tree(text='(A, ((B, C)));')))
  
  expect_equal(c(5L, 5L), NPartitions(list(tree8, ape::unroot(tree8))))
  expect_equal(c(5L, 5L), NPartitions(c(8, 8)))
  expect_error(NPartitions('not a tree'))
})

test_that('Edge distances are calculated correctly', {
  tree <- ape::read.tree(text = '(((a, b), (c, d)), (e, (f, (g, (h, i)))));')
  ed <- EdgeDistances(tree)
  
  expect_equal(ed, t(ed)) # Symmetry
  expect_equal(c(4, 5, 6, 6,
                 5, 6, 6, 4,
                 4, 3, 3, 2, 
                 2, 1, 1, 0),
               ed[, 16])
  expect_equal(c(1, 2, 3, 3, 2, 3, 3, 1, 1, 0, 1, 1, 2, 2, 3, 3), ed[, 10])
  expect_equal(c(1, 0, 1, 1, 1, 2, 2, 1, 2, 2, 3, 3, 4, 4, 5, 5), ed[, 2])
})
