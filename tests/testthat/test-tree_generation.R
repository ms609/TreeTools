context('tree_generation.R')

test_that('Pectinate trees are generated', {
  expect_equal(ape::read.tree(text = '(t1, (t2, (t3, t4)));'),
               PectinateTree(4L))
  expect_equal(ape::read.tree(text = '(a, (b, (c, (d, e))));'),
               PectinateTree(letters[1:5]))
  expect_equal(ape::read.tree(text = '(a, (b, (c, (d, e))));'),
               PectinateTree(ape::read.tree(text = '(a, ((b, c), (d, e)));')))
  expect_equal(ape::read.tree(text = '(Cricocosmia, (Aysheaia, Siberion));'),
               PectinateTree(Lobo.phy[2:4]))
})

test_that('Balanced trees are generated correctly', {
  expect_equal(ape::read.tree(text = '((((t1, t2), t3), (t4, t5)), ((t6, t7), (t8, t9)));'),
               BalancedTree(9L))
  expect_equal(BalancedTree(as.character(1:9)), BalancedTree(1:9))
  escapees <- c("Apostrophe's", 'and quote"s')
  expect_equivalent(PectinateTree(escapees), BalancedTree(escapees))
  expect_equal(integer(0), BalancedBit(seq_len(0)))
  expect_equal('Test', BalancedBit('Test'))
})

test_that("Random trees are generated correctly", {
  expect_equal(c(4, 5, 5, 4, 5), RandomTree(3, root = TRUE)$edge[1:5])
  expect_equal(PectinateTree(c('t2', 't3', 't1')), RandomTree(3, root = 't2'))
  expect_equal(c(4, 4, 4), RandomTree(3, root = FALSE)$edge[1:3])
})

test_that("EnforceOutgroup() fails nicely", {
  expect_error(EnforceOutgroup(BalancedTree(6), 'Non-taxon'))
  expect_equal(BalancedTree(letters[5:6]),
               Subtree(Preorder(EnforceOutgroup(letters[1:8], letters[5:6])), 15))
  expect_equal(ape::root(BalancedTree(8), 't1', resolve.root = TRUE),
               EnforceOutgroup(BalancedTree(8), 't1'))
})
