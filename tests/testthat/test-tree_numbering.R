context("tree_numbering.R")

test_that("RenumberTree handles polytomies", {
  tr <- ape::read.tree(text = '(a, (b, d, c));')
  edge <- tr$edge
  parent <- edge[, 1]
  child <- edge[, 2]

  ret <- RenumberTree(parent, child)
  expect_equal(c(5, 5, 6, 6, 6), ret[, 1])
  expect_equal(c(1, 6, 2, 3, 4), ret[, 2])
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

