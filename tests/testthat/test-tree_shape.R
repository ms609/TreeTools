context('tree_shape.cpp')

test_that('Tree shapes calculated', {
  tree <- PectinateTree(0:3)
  cat("\n\n")
  tree <- Postorder(root(tree, 1, resolve.root = TRUE))
  edge <- tree$edge
  nTip <- NTip(tree)
  plot(tree); nodelabels(4:6); edgelabels(0:5);
  edge <- PostorderEdges(edge[, 1], edge[, 2], nTip = nTip)

  expect_equal(0L, TreeShape(tree))


  PectinateTest <- function (i) expect_equal(0L, TreeShape(PectinateTree(i)))
  lapply(4:16, PectinateTest)
})
