test_that("PaintTree returns correct structure", {
  tree <- BalancedTree(15)
  paint <- PaintTree(tree)

  expect_named(paint, c("edgeCol", "tipCol", "nodeCol"))
  expect_type(paint$edgeCol, "character")
  expect_type(paint$tipCol, "character")
  expect_type(paint$nodeCol, "character")

  expect_length(paint$edgeCol, nrow(tree$edge))
  expect_length(paint$tipCol, NTip(tree))
  expect_length(paint$nodeCol, tree$Nnode)
})

test_that("PaintTree is deterministic", {
  tree <- BalancedTree(8)
  expect_identical(PaintTree(tree), PaintTree(tree))
})

test_that("PaintTree hue partitioning matches specification", {
  tree <- ape::read.tree(text = "(a, (b, (c, d)));")
  # NB this is a binary rooted tree (a, (b, (c, d))) with 4 tips
  # nDesc: a=1, b=1, c=1, d=1; (c,d)=2; (b,(c,d))=3; root=4
  # Root range 0-360; children weights a:1, (b,(c,d)):3
  #   a:        0-90 (mid 45)
  #   (b,(c,d)):90-360 (mid 225)
  # (b,(c,d)) range 90-360; children weights b:1, (c,d):2
  #   b:    90-180 (mid 135)
  #   (c,d):180-360 (mid 270)
  # (c,d) range 180-360; children c:1, d:1
  #   c: 180-270 (mid 225)
  #   d: 270-360 (mid 315)
  paint <- PaintTree(tree, palette = function(h, s) h)
  nTip <- NTip(tree)
  # tipCol is hue at each tip in tip order
  expect_equal(paint$tipCol[match("a", tree$tip.label)], 45)
  expect_equal(paint$tipCol[match("b", tree$tip.label)], 135)
  expect_equal(paint$tipCol[match("c", tree$tip.label)], 225)
  expect_equal(paint$tipCol[match("d", tree$tip.label)], 315)
  # Root is the first internal node (5); (b,(c,d)) is 6; (c,d) is 7
  expect_equal(paint$nodeCol[1], 180)   # root midpoint
  expect_equal(paint$nodeCol[2], 225)   # (b,(c,d)) midpoint
  expect_equal(paint$nodeCol[3], 270)   # (c,d) midpoint
})

test_that("PaintTree saturation invariants", {
  tree <- BalancedTree(8)
  sat <- PaintTree(tree, palette = function(h, s) s)

  # Tips fully saturated
  expect_equal(sat$tipCol, rep_len(1, NTip(tree)))
  # Root achromatic
  expect_equal(sat$nodeCol[1], 0)
  # Saturation strictly increases from root to tip along every edge:
  # for each edge, child's sat >= parent's sat
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  allSat <- c(rep_len(1, NTip(tree)), sat$nodeCol)
  expect_true(all(allSat[child] >= allSat[parent] - 1e-9))
})

test_that("PaintTree root colour is grey under default palette", {
  paint <- PaintTree(BalancedTree(8))
  root <- as.integer(col2rgb(paint$nodeCol[1]))
  expect_equal(root[1], root[2])
  expect_equal(root[2], root[3])
})

test_that("edgeCol always reports the colour of its child", {
  tree <- BalancedTree(10)
  paint <- PaintTree(tree)
  child <- tree$edge[, 2]
  nTip <- NTip(tree)
  isTip <- child <= nTip
  childCol <- character(length(child))
  childCol[isTip] <- paint$tipCol[child[isTip]]
  childCol[!isTip] <- paint$nodeCol[child[!isTip] - nTip]
  expect_equal(paint$edgeCol, childCol)
})

test_that("PaintTree preserves edge identity when $edge is not in preorder", {
  # Build a non-preorder tree by reversing the edge matrix without
  # renumbering any nodes.  Each row of edgeCol must still describe its
  # corresponding row of $edge.
  tree <- BalancedTree(10)
  rev_tree <- tree
  rev_tree$edge <- tree$edge[nrow(tree$edge):1, , drop = FALSE]
  attr(rev_tree, "order") <- NULL

  paint <- PaintTree(rev_tree)
  expect_length(paint$edgeCol, nrow(rev_tree$edge))

  child <- rev_tree$edge[, 2]
  nTip <- NTip(rev_tree)
  isTip <- child <= nTip
  childCol <- character(length(child))
  childCol[isTip] <- paint$tipCol[child[isTip]]
  childCol[!isTip] <- paint$nodeCol[child[!isTip] - nTip]
  expect_equal(paint$edgeCol, childCol)
})

test_that("PaintTree palette dispatch", {
  tree <- BalancedTree(8)
  defaultPaint <- PaintTree(tree, "default")
  protPaint <- PaintTree(tree, "Prot")  # pmatch / tolower
  triPaint <- PaintTree(tree, "TRITANOPIA")

  expect_false(identical(defaultPaint$tipCol, protPaint$tipCol))
  expect_false(identical(defaultPaint$tipCol, triPaint$tipCol))

  redPaint <- PaintTree(tree, function(h, s) rep("#ff0000", length(h)))
  expect_true(all(redPaint$edgeCol == "#ff0000"))
  expect_true(all(redPaint$tipCol == "#ff0000"))
})

test_that("PaintTree rejects bad palette arguments", {
  tree <- BalancedTree(4)
  expect_error(PaintTree(tree, "bogus"), "must match")
  expect_error(PaintTree(tree, c("default", "default")), "single string")
  expect_error(PaintTree(tree, 42), "single string")
})

test_that("PaintTree handles single-tip tree", {
  tree <- structure(list(
    edge = matrix(c(2L, 1L), ncol = 2),
    tip.label = "a",
    Nnode = 1L
  ), class = "phylo")
  paint <- PaintTree(tree)
  expect_length(paint$tipCol, 1)
  expect_length(paint$nodeCol, 1)
})
