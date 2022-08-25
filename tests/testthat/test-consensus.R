test_that("Consensus() errors", {
  expect_error(Consensus(trees), regexp = " a list")
  
  bal8 <- BalancedTree(8)
  oneLeaf <- Consensus(list(DropTip(bal8, 1:7))[c(1, 1, 1)])
  expect_equal(class(oneLeaf), "phylo")
  expect_equal(oneLeaf$Nnode, 0)
  expect_identical(Consensus(list(DropTip(bal8, 1:8))[c(1, 1, 1)]),
                   DropTip(bal8, 1:8))
  expect_warning(expect_true(all.equal(
    Consensus(list(PectinateTree(6), PectinateTree(8))),
    PectinateTree(6)
  )), regexp = "Tree sizes")

  halfTree <- CollapseNode(bal8, 10:12)
  expect_equal(Consensus(halfTree), halfTree)
  expect_equal(Consensus(c(halfTree)), halfTree)
  expect_equal(Consensus(list(halfTree)), halfTree)
})

test_that("Consensus()", {

  ApeTest <- function(tr, p = 1) {
    # plot(Consensus(tr))
    # plot(consensus(tr))
    skip_if_not_installed("ape", "5.5.1") # Bug in ape::consensus?
    if (!expect_true(all.equal(RootTree(Consensus(tr, p = p), 1),
                               RootTree(ape::consensus(tr, p = p), 1)))) {
      dput(RootTree(Consensus(tr, p = p), 1)$edge)
      dput(RootTree(ape::consensus(tr, p = p), 1)$edge)
    }
  }

  trees <- list(BalancedTree(8), PectinateTree(8))[c(1, 1, 1, 1, 2, 2, 2)]
  trees <- RenumberTips(trees, trees[[1]])
  expect_error(Consensus(trees, 2), regexp = "`p`") # p too hign
  expect_error(Consensus(trees, 0.2), regexp = "`p`") # p too low

  ApeTest(trees)
  ApeTest(trees, 0.5)
  expect_equal(Consensus(trees, 0.999), Consensus(trees, 1))

  ApeTest(as.phylo(0:2, 8))
  ApeTest(as.phylo(0:250, 8))
  ApeTest(as.phylo(0:250, 80))
  
  trees <- list(ape::read.tree(text = "((a, b), (c, d));"),
                ape::read.tree(text = "((a, c), (b, d));"))
  expect_equal(Consensus(trees), Preorder(StarTree(letters[1:4])))
})

test_that("ConsensusWithout() is robust", {
  tr <- as.phylo(0:6, 16)
  expect_identical(ConsensusWithout(tr, paste0("t", 1:16)),
                   DropTip(tr[[1]], 1:16))
  expect_identical(ConsensusWithout(tr, paste0("t", 1:4)),
                   Consensus(DropTip(tr, 1:4)))

  expect_true(all.equal(BalancedTree(8), ConsensusWithout(BalancedTree(8))))
  expect_true(all.equal(BalancedTree(4),
                        ConsensusWithout(BalancedTree(8), paste0("t", 5:8))))
  balAndPec <- list(BalancedTree(8), PectinateTree(8))
  t25 <- paste0("t", c(2:5))
  expect_true(all.equal(PectinateTree(paste0("t", c("1", 6:8))),
                        ConsensusWithout(balAndPec, t25)))
  expect_true(all.equal(
    ConsensusWithout(structure(balAndPec, class = "multiPhylo"), t25),
    ConsensusWithout(balAndPec, t25))
  )

  nasty <- structure(list(edge = structure(
    c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
      5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
    .Dim = c(12, 2)),
    Nnode = 5L,
    tip.label = letters[1:8]),
    class = "phylo") # Danger: Do not plot!
  expect_true(all.equal(Preorder(nasty), ConsensusWithout(nasty)))
  expect_true(all.equal(DropTip(nasty, 2), ConsensusWithout(nasty, "b")))

})
