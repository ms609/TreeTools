test_that("SortTree() works", {
  expect_equal(matrix(c(7:10, 10, 9, 11, 11, 8, 7:10, 3:4, 11, 1:2, 5:6), 10),
               SortTree(as.phylo(10, 6))$edge)
  expect_error(#TODO sort unrooted trees,
               SortTree(UnrootTree(PectinateTree(5)))$edge)
  expect_equal('cladewise', attr(SortTree(PectinateTree(5)), 'order'))
})

test_that("SortTree.multiPhylo()", {
  t1 <- as.phylo(123, 12)
  t2 <- as.phylo(921, 12)

  expect_identical(structure(list(SortTree(t1), SortTree(t2)),
                             class = 'multiPhylo'),
                   SortTree(c(t1, t2)))
})
