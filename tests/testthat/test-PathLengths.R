test_that("PathLengths() works", {
  expect_error(PathLengths(matrix(0, 2, 2)),
               "object of class")
  expect_error(PathLengths(BalancedTree(9)),
               "weight")
  bal9 <- BalancedTree(9)
  bal9$edge.length <- 1:16
  pl <- PathLengths(bal9)
  expect_equal(dim(pl), rep(17L, 2))
  Test <- function(tr) {
    calculated <- PathLengths(tr)
    edge <- tr$edge
    parent <- edge[, 1]
    child <- edge[, 2]
    weight <- tr[["edge.length"]]
    nVert <- tr$Nnode + NTip(tr)
    anc <- AllAncestors(parent, child)
    parentEdge <- match(seq_len(nVert), child)
    lapply(seq_len(nVert)[-(NTip(tr) + 1L)], function(end) {
      ancs <- anc[[end]]
      expect_equal(calculated[ancs, end],
                   cumsum(weight[parentEdge[c(end, ancs[-length(ancs)])]]))
    })
  }
  Test(bal9)
  
  pec8 <- CollapseNode(PectinateTree(8), 12:13)
  pec8$edge.length <- 1:12
  Test(pec8)
})
