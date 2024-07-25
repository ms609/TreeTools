test_that("TopologyOnly() works", {
  topol <- as.phylo(10, 5)
  tree <- topol
  tree[["edge.length"]] <- topol[["edge"]][, 1]
  tree[["node.label"]] <- seq_len(topol$Nnode)
  
  # .phylo
  expect_equal(TopologyOnly(tree), topol)
  expect_equal(TopologyOnly(Postorder(tree)), topol)
  
  # .list
  expect_equal(TopologyOnly(list(tree, tree)), list(topol, topol))
  
  # .multiPhylo
  expect_equal(TopologyOnly(c(tree, tree)), c(topol, topol))
  
  # .Splits
  expect_equal(TopologyOnly(as.Splits(tree)), as.Splits(topol))
  
  # .NULL
  expect_null(TopologyOnly(NULL))
})
