test_that("drop_tip() works", {
  expect_error(drop_tip(BalancedTree(8)$edge[, c(1, 2, 1)], 6:4))
  expect_error(drop_tip(BalancedTree(8)$edge[, 1, drop = FALSE], 6:4))
  expect_error(drop_tip(BalancedTree(8)$edge[, 1], 6:4))
  
  expect_equal(matrix(c(6, 7, 8, 8, 7, 6, 9, 9,
                        7, 8, 1, 2, 3, 9, 4, 5), 8, 2),
               drop_tip(BalancedTree(8)$edge, 6:4))
  expect_equal(BalancedTree(4)$edge,
               drop_tip(BalancedTree(8)$edge, 5:8))
  
  testTree <- ape::read.tree(text = "(a, ((b, c), ((d, e, f), g)));")
  testEdge <- Preorder(testTree)$edge
})

test_that("keep_tip() works", {
  expect_error(keep_tip(BalancedTree(8)$edge[, c(1, 2, 1)], !tabulate(6:4, 8)))
  expect_error(keep_tip(BalancedTree(8)$edge[, 1, drop = FALSE],
                        !tabulate(6:4, 8)))
  expect_error(keep_tip(BalancedTree(8)$edge[, 1], !tabulate(6:4, 8)))
  
  expect_equal(keep_tip(BalancedTree(9)$edge, !tabulate(5:8, 9)),
               matrix(c(6, 7, 8, 9, 9, 8, 7, 6,
                        7, 8, 9, 1, 2, 3, 4, 5), 8, 2))
  
  expect_equal(keep_tip(BalancedTree(8)$edge, !tabulate(6:4, 8)),
               matrix(c(6, 7, 8, 8, 7, 6, 9, 9,
                        7, 8, 1, 2, 3, 9, 4, 5), 8, 2))
  expect_equal(keep_tip(BalancedTree(8)$edge, !tabulate(1:4, 8)),
               BalancedTree(4)$edge)
  
  expect_equal(keep_tip(BalancedTree(8)$edge, !tabulate(3:8, 8)),
               BalancedTree(2)$edge)
  
  expect_equal(keep_tip(ape::unroot(BalancedTree(4))$edge, !tabulate(1, 4)),
               matrix(c(4, 4, 4, 1, 2, 3), 3, 2))
  
  unrooted <- ape::read.tree(text = "(a, b, (c, d, ((e1, e2), (f, g))));")
  if(interactive()) {
    plot(unrooted)
    nodelabels()
    tiplabels()
    edgelabels()
  }
  expect_equal(keep_tip(unrooted$edge, !tabulate(1:4, 8)),
               ape::unroot(BalancedTree(4))$edge)
  
  testTree <- ape::read.tree(text = "(a, ((b, c), ((d, e, f), g)));")
  testEdge <- Preorder(testTree)$edge
})

test_that("DropTip() works", {
  bal8 <- BalancedTree(8)
  expect_equal(NTip(DropTip(bal8, 1:8)), 0)
  expect_warning(expect_true(all.equal(bal8, DropTip(bal8, -1))))
  expect_warning(expect_true(all.equal(bal8, DropTip(bal8, 99))))
  expect_warning(expect_true(all.equal(bal8, DropTip(bal8, 'MissingTip'))))
  expect_error(DropTip(bal8, list('Invalid format')))
  
  expect_equal(DropTip(bal8, 7:8), DropTip(bal8, 15L))
  expect_true(all.equal(ape::drop.tip(bal8, 6:8), DropTip(bal8, 6:8)))
  expect_true(all.equal(ape::drop.tip(bal8, c(3, 5, 7)), DropTip(bal8, c(3, 5, 7))))
  
  expect_equal(DropTip(Preorder(nasty), c(1, 3)),
               Preorder(DropTip(nasty, c(1, 3))))
  
  
  bigTree <- RandomTree(1284)
  set.seed(1284)
  bigTip <- sample(1:1284, 608)
  expect_true(all.equal(drop.tip(bigTree, bigTip), DropTip(bigTree, bigTip)))
  #microbenchmark(ape::drop.tip(bigTree, bigTip), DropTip(bigTree, bigTip), times = 25)
  #profvis(replicate(25, DropTip(bigTree, bigTip)), interval = 0.005)
})

test_that("DropTip.multiPhylo() with attributes", {
  multi <- c(bal8 = BalancedTree(8), pec8 = PectinateTree(8))
  attr(multi, 'TipLabel') <- paste0('t', 1:8)
  
  expect_equal(attr(DropTip(multi, 't8'), 'TipLabel'),
               paste0('t', 1:7))
  expect_equal(names(DropTip(multi, 't8')), names(multi))
  expect_equal(DropTip(multi[1], 't1')[[1]], DropTip(multi[[1]], 't1'))
})

test_that("KeepTip() works", {
  expect_warning(expect_true(all.equal(
    BalancedTree(paste0('t', 5:8)),
    KeepTip(BalancedTree(8), paste0('t', 5:9))
  )))
  
  expect_warning(expect_true(all.equal(
    BalancedTree(paste0('t', 5:8)),
    KeepTip(BalancedTree(8), 5:9)
  )))
})
