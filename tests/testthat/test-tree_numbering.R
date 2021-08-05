nastyEdge <- structure(c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
                         5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
                       .Dim = c(12, 2))
nasty <- structure(list(edge = nastyEdge, Nnode = 5L, tip.label = letters[1:8]),
                   class = 'phylo')

test_that("RenumberTree() fails safely", {
  expect_error(RenumberTree(1:3, 1:4))

  Preorder(PectinateTree(8191)) # Largest handled with 16-bit integers
  Preorder(PectinateTree(8192 * 2))

  bigEdge <- PectinateTree(16385)$edge
  expect_error(postorder_edges(bigEdge))
  preorder_edges_and_nodes(bigEdge[, 1], bigEdge[, 2])
})

test_that("RenumberTree() handles polytomies", {
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

  expect_equal(c(9, 10, 10, 11, 11, 11, 10, 12, 12, 13, 13, 9,
                 10, 1, 11,  2,  4,  7, 12,  3, 13,  5,  6, 8),
               as.integer(RenumberTree(nasty$edge[, 1], nasty$edge[, 2])))
})

test_that("RenumberTree() handles singles", {
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

  post6 <- Postorder(BalancedTree(6))$edge
  parent6 <- post6[, 1]
  child6 <- post6[, 2]
  expect_equal(c(9, 9, 11, 11, 8, 8, 10, 10, 7, 7), parent6)
  # Order of tip pairs is arbitrary\
  expect_equal(1:2, sort(child6[parent6 == 9]))
  expect_equal(4:5, sort(child6[parent6 == 11]))
  expect_equal(c(6, 11), sort(child6[parent6 == 10]))
  expect_equal(c(3, 9), sort(child6[parent6 == 8]))
  expect_equal(c(8, 10), sort(child6[parent6 == 7]))

  star <- ape::read.tree(text = '(a, b, d, c);')
  edge <- RenumberTips(star, letters[1:4])$edge
  expect_equal(edge,
               RenumberTips(star, ape::read.tree(text = '(a, b, c, d);'))$edge)
  expect_equal(star$edge, RenumberTree(edge[, 1], edge[, 2]))
  expect_equal(list(star$edge[, 1], star$edge[, 2]),
               RenumberEdges(edge[, 1], edge[, 2]))
})

test_that("RenumberTips() works correctly", {
  abcd <- letters[1:4]
  dcba <- letters[4:1]
  bal7b <- BalancedTree(dcba)
  bal7f <- BalancedTree(abcd)
  pec7f <- PectinateTree(abcd)
  pec7b <- PectinateTree(dcba)

  l7 <- list(bal7b, bal7f, pec7f)
  f7 <- list(bal7f, bal7f, pec7f)
  b7 <- list(bal7b, bal7b, pec7b)
  mp7 <- structure(l7, class = 'multiPhylo')

  expect_equal(f7, RenumberTips(l7, abcd))
  expect_equal(b7, RenumberTips(l7, dcba))

  expect_equal(structure(f7, class = 'multiPhylo'), RenumberTips(mp7, abcd))
  expect_equal(structure(b7, class = 'multiPhylo'), RenumberTips(mp7, dcba))

  expect_error(RenumberTips(l7, letters[1:5]))
  expect_error(RenumberTips(l7, letters[2:5]))
})

test_that("Reorder methods work correctly", {
  bal7 <- BalancedTree(7)
  bal7$edge.length <- 1:12
  attr(bal7, 'order') <- NULL
  pec7 <- PectinateTree(7)
  list7 <- list(bal7, pec7)
  stt <- SingleTaxonTree(1)
  bad <- bal7
  bad$Nnode <- 100
  attr(bad, 'order') <- NULL
  mp7 <- structure(list7, class = 'multiPhylo')
  Test <- function (Method, ..., testEdges = TRUE) {
    expect_identical(Method(bal7, ...), Method(list7, ...)[[1]])
    expect_identical(Method(pec7, ...), Method(mp7, ...)[[2]])
    expect_equal(stt, Method(stt))
    expect_identical(Method(bal7), Method(Method(bal7)))
    if (testEdges) expect_equal(Method(bal7)$edge, Method(bal7$edge))
    expect_error(Method(10))
    expect_error(Method(1:2))
    expect_error(Method(matrix('one')))
  }
  Test(ApePostorder, testEdges = FALSE)
  expect_error(ApePostorder(bad))

  Test(Postorder)
  expect_equal(c(rep(13:11, each = 2), 11, 10, 10, 10, 9, 9),
               Postorder(nastyEdge, renumber = TRUE)[, 1])

  Test(Cladewise)
  expect_error(Cladewise(bad))

  Test(Preorder)

  Test(Pruningwise, testEdges = FALSE)
  expect_error(Pruningwise(bad))

})

test_that("Malformed trees don't cause crashes", {
  treeDoubleNode <- read.tree(text = "((((((1,2)),3),4),5),6);")
  treePolytomy   <- read.tree(text = "((((1,2,3),4),5),6);")
  treeDoublyPoly <- read.tree(text = "(((((1,2,3)),4),5),6);")
  nasty <- structure(list(edge = structure(
    c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
      5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
    .Dim = c(12, 2)),
    Nnode = 5L,
    tip.label = letters[1:8]),
    class = 'phylo') # Danger: Do not plot!

  reordered <- Preorder(treeDoubleNode)$edge
  expect_equal(11L, dim(reordered)[1])
  expect_equal(5L, sum(tabulate(reordered[, 1]) == 2L))

  postordered <- Postorder(treeDoubleNode)$edge
  expect_equal(11L, dim(postordered)[1])
  expect_equal(5L, sum(tabulate(postordered[, 1]) == 2L))


  reordered <- Preorder(treePolytomy)$edge
  expect_equal(9L, dim(reordered)[1])
  expect_equal(c(2L, 2L, 2L, 3L), as.integer(table(reordered[, 1])))

  reordered <- Postorder(treePolytomy)$edge
  expect_equal(9L, dim(reordered)[1])
  expect_equal(c(2L, 2L, 2L, 3L), as.integer(table(reordered[, 1])))


  reordered <- Preorder(treeDoublyPoly)$edge
  expect_equal(10L, dim(reordered)[1])
  expect_equal(c(2L, 2L, 2L, 1L, 3L), as.integer(table(reordered[, 1])))

  reordered <- Postorder(treeDoublyPoly)$edge
  expect_equal(10L, dim(reordered)[1])
  expect_equal(c(2L, 2L, 2L, 1L, 3L), as.integer(table(reordered[, 1])))

  #C <- 0
  #plot(Preorder(nasty)); nodelabels(c(12, 10, 13, 11, 9) - C); tiplabels(1:8 - C)
  #edgelabels(c(2, 3, 6, 4, 8, 10, 9, 5, 7, 1, 12, 11) - C)
  reordered <- Preorder(nasty)$edge
  expect_equal(12L, dim(reordered)[1])
  # Nodes renumbered
  expect_equal(c(2L, 3L, 3L, 2L, 2L), tabulate(reordered[, 1])[9:13])

  reordered <- Postorder(nasty)$edge
  expect_equal(12L, dim(reordered)[1])
  # Original node numbers retained
  expect_equal(c(2L, 3L, 2L, 2L, 3L), tabulate(reordered[, 1])[9:13])

})

