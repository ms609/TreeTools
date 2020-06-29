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
  expect_true(is.integer(PectinateTree(8)$edge))
})

test_that('Balanced trees are generated correctly', {
  # nTip even
  expect_equal(ape::read.tree(text = '(((t1, t2), (t3, t4)), ((t5, t6), (t7, t8)));'),
               BalancedTree(8L))
  # nTip odd
  expect_equal(ape::read.tree(text = '((((t1, t2), t3), (t4, t5)), ((t6, t7), (t8, t9)));'),
               BalancedTree(9L))
  expect_equal(BalancedTree(as.character(1:9)), BalancedTree(1:9))
  escapees <- c("Apostrophe's", 'and quote"s')
  expect_equivalent(PectinateTree(escapees), BalancedTree(escapees))
  expect_equal(integer(0), .BalancedBit(seq_len(0)))
  expect_equal('Test', .BalancedBit('Test'))
  expect_true(is.integer(BalancedTree(8)$edge))
})

test_that("StarTree() works", {
  expect_equal(ape::read.tree(text = '(t1, t2, t3, t4, t5, t6, t7, t8);'),
               StarTree(8L))
  expect_true(is.integer(StarTree(8)$edge))
})

test_that("Random trees are generated correctly", {
  expect_equal(c(4, 4, 5, 5, 1, 5, 2, 3), RandomTree(3, root = TRUE)$edge[1:8])
  expect_equal(PectinateTree(c('t2', 't3', 't1')), RandomTree(3, root = 't2'))
  expect_equal(c(4, 4, 4), RandomTree(3, root = FALSE)$edge[1:3])
})

test_that("NJTree() works", {
  a..f <- letters[1:6]
  bal6 <- StringToPhyDat('111100 111000 111000 110000', letters[1:6],
                         byTaxon = FALSE)
  expect_equal(BalancedTree(letters[c(1:3, 6:4)]), RootTree(NJTree(bal6), a..f[1:3]))
  expect_equal(c(0, 1, 2, 1, rep(0, 6)),
               Preorder(NJTree(bal6, TRUE))$edge.length * 4L)
})

test_that("EnforceOutgroup() fails nicely", {
  expect_error(EnforceOutgroup(BalancedTree(6), 'Non-taxon'))
  expect_error(EnforceOutgroup(BalancedTree(6), c('t1', 'Non-taxon')))
  expect_equal(BalancedTree(letters[5:6]),
               Subtree(Preorder(EnforceOutgroup(letters[1:8], letters[5:6])), 15))
  expect_equal(ape::root(BalancedTree(8), 't1', resolve.root = TRUE),
               EnforceOutgroup(BalancedTree(8), 't1'))
})

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
  expect_error(.RandomParent(1))
  expect_error(.RandomParent(0))
  expect_error(.RandomParent(-1))

  nSamples <- 100
  # Ape's trees are not uniformly distributed:
  expect_lt(chisq.test(table(vapply(lapply(rep(5, nSamples), ape::rtree),
                    as.TreeNumber, 1)))$p.value,
            0.001)

  # Our trees are:
  expect_gt(chisq.test(table(vapply(lapply(rep(5, nSamples), RandomTree),
                    as.TreeNumber, 1)))$p.value,
            0.001)

  expect_false(is.rooted(RandomTree(6, root = FALSE)))
  expect_true(is.rooted(RandomTree(6, root = TRUE)))
})
