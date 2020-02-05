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

test_that("RenumberTree handles singles", {
  withSingles <- ape::read.tree(text='(a, (b, (c), (((d), (e)))));')
  expect_equal(c(6, 6, 7, 7, 8, 7, 9, 10, 11, 10, 12,
                 1, 7, 2, 8, 3, 9, 10, 11, 4, 12, 5),
               as.integer(Preorder(withSingles)$edge))
})

test_that("Replacement reorder functions work correctly", {
  ## Tree
  tree <- ape::read.tree(text = "((((((1,2),3),4),5),6),(7,(8,(9,(10,(11,12))))));")
  expect_equal(ape::reorder.phylo(tree, 'cladewise'), Cladewise(tree))
  expect_equal(ape::reorder.phylo(tree, 'pruningwise'), Pruningwise(tree))

  expect_equal(matrix(c(9,9,11,11,10,10,8,8,7,7,
                        2,3,5,6,4,11,1,9,8,10), ncol=2),
               Postorder(BalancedTree(6))$edge)


  star <- ape::read.tree(text = '(a, b, d, c);')
  edge <- RenumberTips(star, letters[1:4])$edge
  expect_equal(edge,
               RenumberTips(star, ape::read.tree(text = '(a, b, c, d);'))$edge)
  expect_equal(star$edge, RenumberTree(edge[, 1], edge[, 2]))
  expect_equal(list(star$edge[, 1], star$edge[, 2]),
               RenumberEdges(edge[, 1], edge[, 2]))
})

test_that("Preorder handles malformed trees without crashing", {
  treeDoubleNode <- read.tree(text = "((((((1,2)),3),4),5),6);")
  treePolytomy   <- read.tree(text = "((((1,2,3),4),5),6);")
  treeDoublyPoly <- read.tree(text = "(((((1,2,3)),4),5),6);")
  
  reordered <- Preorder(treeDoubleNode)$edge
  expect_equal(11L, dim(reordered)[1])
  expect_equal(5L, sum(table(reordered[, 1]) == 2L))
  
  reordered <- Preorder(treePolytomy)$edge
  expect_equal(9L, dim(reordered)[1])
  expect_equal(c(2L, 2L, 2L, 3L), as.integer(table(reordered[, 1])))
  
  reordered <- Preorder(treeDoublyPoly)$edge
  expect_equal(10L, dim(reordered)[1])
  expect_equal(c(2L, 2L, 2L, 1L, 3L), as.integer(table(reordered[, 1])))
})