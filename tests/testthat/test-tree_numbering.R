context("tree_numbering.R")

test_that("RenumberTree handles polytomies", {
  tr <- ape::read.tree(text = '(a, (b, d, c));')
  edge <- tr$edge
  parent <- edge[, 1]
  child <- edge[, 2]

  ret <- RenumberTree(parent, child)
  expect_equal(c(5, 5, 6, 6, 6), ret[, 1])
  expect_equal(c(1, 6, 2, 3, 4), ret[, 2])


  edge <- structure(c(6L, 7L, 5L, 7L, 6L, 5L,
                      2L, 5L, 3L, 6L, 1L, 4L),
                    .Dim = c(6L, 2L))

  # Must be in preorder; i.e. each node in left subtree before each node in
  # right subtree for each subtree
  #
  # Also, nodes should be rotated such that the lowest tip in a subtree
  # is always encountered first.
  #
  # These rules ensure a unique representation for any tree.
  expectation <- structure(c(5L, 6L, 6L, 5L, 7L, 7L,
                             6L, 1L, 2L, 7L, 3L, 4L),
                           .Dim = c(6L, 2L))
  expect_equal(expectation,
               RenumberTree(edge[, 1], edge[, 2]))

  nasty <- structure(list(edge = structure(
    c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
      5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
                                 .Dim = c(12, 2)),
                Nnode = 5L,
                tip.label = letters[1:8]),
                class = 'phylo')
  expect_equal(c(9, 10, 10, 11, 11, 11, 10, 12, 12, 13, 13, 9,
                 10, 1, 11,  2,  4,  7, 12,  3, 13,  5,  6, 8),
               as.integer(RenumberTree(nasty$edge[, 1], nasty$edge[, 2])))

})

test_that("replacement reorder functions work correctly", {
  ## Tree
  tree <- ape::read.tree(text = "((((((1,2),3),4),5),6),(7,(8,(9,(10,(11,12))))));")
  expect_equal(ape::reorder.phylo(tree, 'cladewise'), Cladewise(tree))
  expect_equal(ape::reorder.phylo(tree, 'pruningwise'), Pruningwise(tree))
  expect_equal(ape::reorder.phylo(tree, 'postorder'), Postorder(tree))

  star <- ape::read.tree(text = '(a, b, d, c);')
  edge <- RenumberTips(star, letters[1:4])$edge
  expect_equal(star$edge, RenumberTree(edge[, 1], edge[, 2]))
  expect_equal(list(star$edge[, 1], star$edge[, 2]),
               RenumberEdges(edge[, 1], edge[, 2]))
})

