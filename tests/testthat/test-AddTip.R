nasty <- structure(list(edge = structure(
  c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
    5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
  .Dim = c(12, 2)),
  Nnode = 5L,
  tip.label = letters[1:8]),
  class = "phylo") # Danger: Do not plot!

test_that("AddTipEverywhere() correct", {
  backbone <- PectinateTree(5)
  
  expect_equal(7L, length(AddTipEverywhere(backbone, includeRoot = FALSE)))
  expect_equal(7L, length(AddTipEverywhere(UnrootTree(backbone),
                                           includeRoot = FALSE)))
  expect_equal(9L, length(AddTipEverywhere(backbone, includeRoot = TRUE)))
  expect_equal(5L, length(AddTipEverywhere(CollapseNode(backbone, 7:9))))
  
})

test_that("AddTip() at root", {
  expect_true(all.equal(
    AddTip(DropTip(PectinateTree(8), 1), where = 0, label = "t1"),
    PectinateTree(8)))
  AddTip(nasty, 1L)
  AddTip(nasty, 12L)
  expect_true(TRUE) # no crashes!
})

test_that("AddTip() with tip name", {
  bal8 <- BalancedTree(8)
  expect_error(AddTip(bal8, "invalid tip"))
  expect_equal(AddTip(bal8, 1L), AddTip(bal8, "t1"))
  expect_equal(AddTip(nasty, 1L), AddTip(nasty, "a"))
})

test_that("AddTip() with edge lengths", {
  pec8 <- PectinateTree(8)
  pec8$edge.length <- rep(1L, 14L)
  
  add1 <- AddTip(pec8, 0, "Root")
  expect_equal(dim(add1$edge)[1], length(add1[["edge.length"]]))
  
  # case = 1 -> y is bound on the root of x
  expect_equal(
    AddTip(pec8, 0, edgeLength = NULL, lengthBelow = 0.6)$edge.length,
    c(0.6, rep(1, 14), 0.6)
  )
  expect_equal(
    AddTip(pec8, 0, edgeLength = 2, lengthBelow = 0.6)$edge.length,
    c(0.6, rep(1, 14), 2)
  )
  # case = 2 -> y is bound on a tip of x
  expect_equal(
    AddTip(pec8, 2, edgeLength = NULL, lengthBelow = 0.7)$edge.length,
    c(rep(1, 2), 0.3, 0.7, 0.7, rep(1, 11))
  )
  expect_equal(
    AddTip(pec8, 2, edgeLength = 2, lengthBelow = 0.7)$edge.length,
    c(rep(1, 2), 0.3, 0.7, 2, rep(1, 11))
  )
  expect_equal(
    AddTip(pec8, 2, edgeLength = NULL, lengthBelow = NULL)$edge.length,
    c(rep(1, 2), 0.5, 0.5, 0.5, rep(1, 11))
  )
  expect_equal(
    AddTip(pec8, 2, edgeLength = 2, lengthBelow = NULL)$edge.length,
    c(rep(1, 2), 0.5, 0.5, 2, rep(1, 11))
  )
  # case = 3 -> y is bound on a node of x
  expect_equal(
    AddTip(pec8, 15, edgeLength = NULL, lengthBelow = 0.7)$edge.length,
    c(rep(1, 11), 0.3, 0.7, 0.7, 1, 1)
  )
  expect_equal(
    AddTip(pec8, 15, edgeLength = 2, lengthBelow = 0.7)$edge.length,
    c(rep(1, 11), 0.3, 2, 0.7, 1, 1)
  )
  expect_equal(
    AddTip(pec8, 15)$edge.length,
    c(rep(1, 11), 0.5, 0, 0.5, 1, 1)
  )
  expect_equal(
    AddTip(pec8, 15, edgeLength = NULL, lengthBelow = 1)$edge.length,
    c(rep(1, 11), 0, 1, 1, 1, 1)
  )
  expect_equal(
    AddTip(pec8, 15, edgeLength = 2, lengthBelow = 1)$edge.length,
    c(rep(1, 11), 0, 2, 1, 1, 1)
  )
  expect_equal(
    AddTip(pec8, 15, edgeLength = NULL, lengthBelow = 1.6)$edge.length,
    c(rep(1, 11), -0.6, 1.6, 1.6, 1, 1)
  )
  expect_equal(
    AddTip(pec8, 15, edgeLength = 2, lengthBelow = 1.6)$edge.length,
    c(rep(1, 11), -0.6, 2, 1.6, 1, 1)
  )
})

test_that("AddTip() handles node labels", {
  bal6 <- BalancedTree(6)
  startLabels <- paste("Node", 7:11)
  bal6[["node.label"]] <- startLabels
  if (interactive()) {
    plot(bal6, show.node.label = TRUE)
  }
  
  # case = 1 -> y is bound on the root of x
  expect_equal(
    AddTip(bal6, where = 0)[["node.label"]],
    c("", startLabels)
  )
  
  # case = 2 -> y is bound on a tip of x
  expect_equal(
    AddTip(bal6, where = 1)[["node.label"]],
    c(startLabels[1:3], "", startLabels[4:5])
  )
  
  # case = 3 -> y is bound on a node of x
  expect_equal(
    AddTip(bal6, where = 11)[["node.label"]],
    c(startLabels[1:4], "", startLabels[[5]])
  )
  
  expect_equal(AddTipEverywhere(bal6)[[1]][["node.label"]],
               AddTip(bal6, where = 1)[["node.label"]])
})

test_that("AddTip(lengthBelow = NA)", {
  tree <- BalancedTree(10)
  tree$edge.length <- 1 + (1:18 / 100)
  tree$node.label <- paste("n", 11:19)
  
  # Case 1: At root
  at11 <- AddTip(tree, 11, "NEW_TIP", lengthBelow = NA)
  expect_equal(at11$edge, rbind(tree$edge + ifelse(tree$edge > 10, 1, 0),
                                c(12, 11)))
  expect_equal(at11$edge.length, c(tree$edge.length, 0))
  expect_equal(at11$node.label, tree$node.label)
  
  # Case 2: At leaf
  at5 <- AddTip(tree, 5, "NEW_TIP", lengthBelow = NA)
  new5 <- tree$edge + ifelse(tree$edge > 10, 1, 0)
  expect_equal(at5$edge, AddTip(tree, 5)$edge)
  expect_equal(at5$edge.length[-10:-11], tree$edge.length)
  expect_equal(at5$edge.length[10:11], c(0, 0))
  expect_equal(at5$node.label[-6], tree$node.label)
  expect_equal(at5$node.label[6], "")
  
  # Case 3: Internal node
  at15 <- AddTip(tree, 15, "NEW_TIP", lengthBelow = NA)
  expect_equal(at15$edge[-8, ], tree$edge + ifelse(tree$edge > 10, 1, 0))
  expect_equal(at15$edge[8, ], c(16, 11))
  expect_equal(at15$edge.length[-8], tree$edge.length)
  expect_equal(at15$edge.length[8], 0)
  expect_equal(at15$node.label, tree$node.label)
})

test_that("AddTipEverywhere() handles nasty tree", {
  added <- AddTipEverywhere(nasty)
  lapply(added, function(tr) expect_true(all(tr$edge > 0)))
  expect_true(all.equal(lapply(added, Preorder),
                        AddTipEverywhere(Preorder(nasty))))
})

test_that("AddTipEverywhere() with tiny trees", {
  added <- AddTipEverywhere(StarTree(2))
  lapply(added, function(tr) expect_true(all(tr$edge > 0)))
  expect_equal(2, length(added))
  expect_equal(3, length(AddTipEverywhere(StarTree(2), include = TRUE)))
  
  expect_equal(list(PectinateTree(c("t1", "New tip"))),
               AddTipEverywhere(StarTree(1)))
  expect_equal(list(SingleTaxonTree("New tip")),
               AddTipEverywhere(structure(list(tip.label = character(0)),
                                          class = "phylo")))
})
