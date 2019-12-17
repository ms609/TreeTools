context('tree_shape.cpp')

test_that('Tree shapes counted', {
  expect_equal(c(1, 1, 1, 2, 3, 6, 11, 23),
               vapply(1:8, NRootedShapes, 0L))

  expect_equal(c(1, 1, 1, 1, 1, 2, 2, 4, 6),
               vapply(1:9, NUnrootedShapes, 0L))
})

test_that('Rooted tree shapes calculated', {
  expect_equal(NRootedShapes(8) - 1L, RootedTreeShape(BalancedTree(0:7)))

  expect_equal(0L, RootedTreeShape(PectinateTree(0:3)))
  expect_equal(0L, UnrootedTreeShape(PectinateTree(0:3)))
  expect_equal(1L, RootedTreeShape(BalancedTree(0:3)))

  expect_equal(0L, RootedTreeShape(PectinateTree(0:4)))
  expect_equal(0L, UnrootedTreeShape(PectinateTree(0:4)))
  expect_equal(1L, RootedTreeShape(as.phylo(1, 5)))
  expect_equal(2L, RootedTreeShape(BalancedTree(0:4)))

  expect_equal(0L, RootedTreeShape(PectinateTree(0:5)))
  expect_equal(0L, UnrootedTreeShape(PectinateTree(0:5)))
  expect_equal(1L, RootedTreeShape(as.phylo(58, 6)))
  expect_equal(2L, RootedTreeShape(as.phylo(1, 6)))
  expect_equal(3L, RootedTreeShape(ape::read.tree(text='((a, b), (c, (d, (e, f))));')))
  expect_equal(4L, RootedTreeShape(ape::read.tree(text='((a, b), ((c, d), (e, f)));')))
  expect_equal(5L, RootedTreeShape(BalancedTree(0:5)))

  PectinateTest <- function (i) expect_equal(0L, RootedTreeShape(PectinateTree(i)))
  lapply(4:16, PectinateTest)

  BalancedTest <- function (i) expect_equal(NRootedShapes(i) - 1L,
                                            RootedTreeShape(BalancedTree(i)))
  lapply(c(2^(1:4), 10), BalancedTest)
})


test_that('Rooted tree shapes built', {
  expect_equal(RootedTreeWithShape(0, 4), PectinateTree(rep('', 4)))
  expect_equal(RootedTreeWithShape(1, 4), BalancedTree(rep('', 4)))

  expect_equal(RootedTreeWithShape(0L, 5L), PectinateTree(rep('', 5)))
  expect_equal(RootedTreeWithShape(1L, 5L), as.phylo(1, 5, rep('', 5)))
  expect_equal(RootedTreeWithShape(2L, 5L), BalancedTree(rep('', 5)))

  blank6 <- rep('', 6L)
  expect_equal(RootedTreeWithShape(0L, 6L), PectinateTree(blank6))
  expect_equal(RootedTreeWithShape(1L, 6L), as.phylo(58, 6, blank6))
  expect_equal(RootedTreeWithShape(2L, 6L), as.phylo(1, 6, blank6))
  expect_equal(RootedTreeWithShape(3L, 6L), ape::read.tree(text='((,),(,(,(,))));'))
  expect_equal(RootedTreeWithShape(4L, 6L), ape::read.tree(text='((,),((,),(,)));'))
  expect_equal(RootedTreeWithShape(5L, 6L), BalancedTree(blank6))

  expect_equal(RootedTreeWithShape(0, 8), PectinateTree(rep('', 8)))
  expect_equal(RootedTreeWithShape(NRootedShapes(8), 8), BalancedTree(rep('', 8)))

  BalancedTest <- function (i) {
    expect_equal(BalancedTree(rep('', i)), RootedTreeWithShape(NRootedShapes(i) - 1L, i))
  }
  lapply(2^(1:4), BalancedTest)
})

test_that('Unrooted tree shapes built', {
  expect_error(plot(UnrootedTreeWithShape(4, 8))) # Out of range

  expect_equal(UnrootedTreeWithShape(0, 9), unroot(PectinateTree(rep('', 9))))
  TestSym <- function (tree, shape) {
    expect_equal(shape, UnrootedTreeShape(tree))
    expect_equal(UnrootedTreeKey(tree),
                 UnrootedTreeKey(UnrootedTreeWithShape(shape, NTip(tree))))
  }
  blank9 <- rep('', 9)
  TestSym(PectinateTree(blank9), 0)
  TestSym(as.phylo(72292, nTip = 9, blank9), 1)
  TestSym(as.phylo(67987, nTip = 9, blank9), 2)
  TestSym(BalancedTree(blank9), 3)
  TestSym(as.phylo(72237, nTip = 9, blank9), 3)
  TestSym(as.phylo(67882, nTip = 9, blank9), 4)
  TestSym(as.phylo(72298, nTip = 9, blank9), 5)

})
