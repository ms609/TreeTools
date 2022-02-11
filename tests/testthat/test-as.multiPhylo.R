test_that("as.multiPhylo()", {
  expect_equal(structure(list(BalancedTree(8), PectinateTree(8)),
                         class = 'multiPhylo'),
               as.multiPhylo(list(BalancedTree(8), PectinateTree(8))))
  expect_equal(structure(list(BalancedTree(8)), class = 'multiPhylo'),
               as.multiPhylo(BalancedTree(8)))

  char <- MatrixToPhyDat(matrix(c(1,1,1,0,0,0), ncol = 1,
                                dimnames = list(letters[1:6], NULL)))
  expect_equal(ape::read.tree(text = '((a, b, c), (d, e, f));'),
               as.multiPhylo(char)[[1]])

  char2 <- MatrixToPhyDat(matrix(c(1,1,1,0,0,0,
                                   0,1,1,1,1,'?',
                                   0,0,1,1,'-',2,
                                   0,1,1,1,1,'?'), ncol = 4,
                                dimnames = list(letters[1:6], NULL)))
  mpChar2 <- as.multiPhylo(char2)

  expect_equal(ape::read.tree(text = '((a, b, c), (d, e, f));'),
               mpChar2[[1]])
  expect_equal(ape::read.tree(text = '(b, c, d, e);'),
               mpChar2[[2]])
  expect_equal(ape::read.tree(text = '((a, b), (c, d));'),
               mpChar2[[3]])
  expect_equal(mpChar2[[2]], mpChar2[[4]])

  mpSplits <- as.Splits(PectinateTree(letters[1:6]))
  expect_true(all.equal(as.multiPhylo(mpSplits), structure(list(
    '9' = ape::read.tree(text = '((a, b), (c, d, e, f));'),
    '10' = ape::read.tree(text = '((a, b, c), (d, e, f));'),
    '11' = ape::read.tree(text = '((a, b, c, d), (e, f));')),
    class = 'multiPhylo')))

})
