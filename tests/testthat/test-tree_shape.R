test_that("Errors are handled", {
  skip_if(Sys.getenv("USING_ASAN") != "")
  expect_error(RootedTreeWithShape(as.integer64(-1)),
               "Shape may not be negative")
  expect_error(RootedTreeWithShape(as.integer64(0), -1),
               "Tree must have at least zero tips")
  expect_error(RootedTreeWithShape(as.integer64(2)^31),
               "Shapes this large are not.* implemented")
  expect_error(UnrootedTreeWithShape(31, 31), " < 31 leaves")
  expect_error(.UnrootedKeys(29L), " 29 leaves")

    expect_error(NRootedShapes(56L), "Too many shapes ")
  expect_error(NUnrootedShapes(61L), "Too many shapes ")
})

test_that("Tree shapes counted", {
  expect_equal(as.integer64(c(1, 1, 1, 2, 3, 6, 11, 23)),
               structure(vapply(1:8, NRootedShapes, bit64::integer64(1L)),
                         class = "integer64"))

  expect_equal(as.integer64(lengths(unrootedKeys)),
               structure(vapply(seq_along(unrootedKeys),
                                NUnrootedShapes,
                                integer64(1)),
                         class = "integer64"))
})

test_that("Nasty node order not fatal", {
  nastyBinary <- structure(list(edge = structure(
    c(9, 12, 10, 13, 11, 10, 11, 14, 15, 13, 12, 9, 14, 15,
      5, 10, 15, 14,  3, 13,  9,  4, 11,  7,  8, 6,  2,  1),
    .Dim = c(14, 2)),
    Nnode = 7L,
    tip.label = letters[1:8]),
    class = "phylo") # Danger: Do not plot!
  expect_equal(RootedTreeShape(Preorder(nastyBinary)),
               RootedTreeShape(nastyBinary))

  skip_if(Sys.getenv("USING_ASAN") != "")
  expect_error(edge_to_rooted_shape(1:10, 1:11, 6), "must be the same length")
  expect_error(edge_to_rooted_shape(1:10, 1:10, 5), "is tree binary\\?")
})

test_that("Rooted tree shapes counted", {
  # https://oeis.org/A001190
  expect_equal(as.integer64(c(1, 1, 1, 2, 3, 6, 11, 23, 46, 98, 207, 451, 983,
                              2179, 4850, 10905, 24631, 56011, 127912, 293547,
                              676157, 1563372, 3626149, 8436379, 19680277,
                              46026618, 107890609, 253450711, 596572387,
                              1406818759, "3323236238",
                              "7862958391", "18632325319", "44214569100")),
               structure(vapply(1:34, NRootedShapes, integer64(1)),
                                           class = "integer64"))
})

test_that("Rooted tree shapes fail gracefully", {
  skip_if(Sys.getenv("USING_ASAN") != "")
  expect_error(RootedTreeShape(BalancedTree(56)), "> 55 leaves")
})

test_that("Rooted tree shapes calculated", {
  expect_equal(NRootedShapes(8) - 1L, RootedTreeShape(BalancedTree(0:7)))
  expect_equal(NRootedShapes(54) - 1L, RootedTreeShape(BalancedTree(54)))
  expect_equal(NRootedShapes(55), RootedTreeShape(BalancedTree(55)))

  expect_equal(as.integer64(0L), RootedTreeShape(PectinateTree(0:3)))
  expect_equal(as.integer64(0L), RootedTreeShape(SortTree(PectinateTree(0:3))))
  expect_equal(0L, UnrootedTreeShape(PectinateTree(0:3)))
  expect_equal(as.integer64(1L), RootedTreeShape(BalancedTree(0:3)))

  expect_equal(as.integer64(0L), RootedTreeShape(PectinateTree(0:4)))
  expect_equal(0L, UnrootedTreeShape(PectinateTree(0:4)))
  expect_equal(as.integer64(1L), RootedTreeShape(as.phylo(1, 5)))
  expect_equal(as.integer64(2L), RootedTreeShape(BalancedTree(0:4)))

  expect_equal(as.integer64(0L), RootedTreeShape(PectinateTree(0:5)))
  expect_equal(0L, UnrootedTreeShape(PectinateTree(0:5)))
  expect_equal(as.integer64(1L), RootedTreeShape(as.phylo(58, 6)))
  expect_equal(as.integer64(2L), RootedTreeShape(as.phylo(1, 6)))
  expect_equal(as.integer64(3L), RootedTreeShape(
    ape::read.tree(text = "((a, b), (c, (d, (e, f))));")))
  expect_equal(as.integer64(4L), RootedTreeShape(
    ape::read.tree(text = "((a, b), ((c, d), (e, f)));")))
  expect_equal(as.integer64(5L), RootedTreeShape(BalancedTree(0:5)))

  PectinateTest <- function(i) {
    expect_equal(as.integer64(0L), RootedTreeShape(PectinateTree(i)))
  }
  lapply(4:16, PectinateTest)

  BalancedTest <- function(i) {
    expect_equal(NRootedShapes(i) - 1L,
                 RootedTreeShape(BalancedTree(i)))
  }
  lapply(c(2 ^ (1:4), 10), BalancedTest)

  expect_equal(0L, .UnrootedKeys(4))
  expect_equal(0L, .UnrootedKeys(5))
  expect_equal(0:1, .UnrootedKeys(6))
  expect_equal(0:1, .UnrootedKeys(7))
  expect_equal(c(0:2, 4), .UnrootedKeys(8))
  expect_equal(c(0:2, 4:5, 7), .UnrootedKeys(9))
})

expect_treequal <- function(...) expect_true(all.equal(...))
test_that("Rooted tree shapes built", {
  expect_error(RootedTreeWithShape(-1, 5))

  expect_treequal(RootedTreeWithShape(0, 4), PectinateTree(rep("", 4)))
  expect_treequal(RootedTreeWithShape(1, 4), BalancedTree(rep("", 4)))

  expect_treequal(RootedTreeWithShape(0L, 5L), PectinateTree(rep("", 5)))
  expect_treequal(RootedTreeWithShape(1L, 5L), as.phylo(1, 5, rep("", 5)))
  expect_treequal(RootedTreeWithShape(2L, 5L), BalancedTree(rep("", 5)))

  blank6 <- rep("", 6L)
  expect_treequal(RootedTreeWithShape(0L, 6L), PectinateTree(blank6))
  expect_treequal(RootedTreeWithShape(1L, 6L), as.phylo(58, 6, blank6))
  expect_treequal(RootedTreeWithShape(2L, 6L), as.phylo(1, 6, blank6))
  expect_treequal(RootedTreeWithShape(3L, 6L),
                  ape::read.tree(text = "((,),(,(,(,))));"))
  expect_treequal(RootedTreeWithShape(4L, 6L),
                  ape::read.tree(text = "((,),((,),(,)));"))
  expect_treequal(RootedTreeWithShape(5L, 6L), BalancedTree(blank6))

  expect_treequal(RootedTreeWithShape(0, 8), PectinateTree(rep("", 8)))
  expect_treequal(RootedTreeWithShape(NRootedShapes(8), 8),
                  BalancedTree(rep("", 8)))

  BalancedTest <- function(i) {
    expect_treequal(BalancedTree(rep("", i)),
                    RootedTreeWithShape(NRootedShapes(i) - 1L, i))
  }
  lapply(2 ^ (1:4), BalancedTest)
})

test_that("Unrooted tree shapes fail gracefully", {
  skip_if(Sys.getenv("USING_ASAN") != "")
  expect_error(UnrootedTreeWithShape(4, 8), "must be between 0 and 3")
})

test_that("Unrooted tree shapes built", {
  expect_treequal(UnrootedTreeWithShape(0, 9),
                  UnrootTree(PectinateTree(rep("", 9))))
  TestSym <- function(tree, shape) {
    expect_equal(shape, UnrootedTreeShape(tree))
    expect_equal(UnrootedTreeKey(tree),
                 UnrootedTreeKey(UnrootedTreeWithShape(shape, NTip(tree))))
  }
  blank9 <- rep("", 9)
  TestSym(PectinateTree(blank9), 0)
  TestSym(as.phylo(72292, nTip = 9, blank9), 1)
  TestSym(as.phylo(67987, nTip = 9, blank9), 2)
  TestSym(BalancedTree(blank9), 3)
  TestSym(as.phylo(72237, nTip = 9, blank9), 3)
  TestSym(as.phylo(67882, nTip = 9, blank9), 4)
  TestSym(as.phylo(72298, nTip = 9, blank9), 5)

})
