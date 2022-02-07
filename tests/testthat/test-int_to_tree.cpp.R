test_that("Failures are graceful", {
  expect_error(num_to_parent(10, 1))
  expect_error(num_to_parent(10, -1))
  expect_error(mixed_base_to_parent(10, 1))
  expect_error(mixed_base_to_parent(10, -1))
  expect_error(edge_to_num(1:10, 1:11, 6))
  expect_error(edge_to_num(1:10, 1:10, 5))
  expect_error(edge_to_mixed_base(1:10, 1:11, 6))
  expect_error(edge_to_mixed_base(1:10, 1:10, 5))
  expect_error(as.phylo(0, 0))
})

test_that("Edge cases handled", {
  expect_equal(SingleTaxonTree('singleton'), as.phylo(0, 1, 'singleton'))
  expect_true(all.equal(BalancedTree(2), as.phylo(0, 2)))
  expect_equal(structure(as.integer64(0L), nTip = 2L, tip.label = c('t1', 't2'),
                         class = c('TreeNumber', 'integer64')),
               as.TreeNumber(as.phylo(0, 2)))
  expect_equal(structure(integer64(1), nTip = 1L, tip.label = 't1',
                         class = c('TreeNumber', 'integer64')),
               as.TreeNumber(SingleTaxonTree()))
})

test_that('Trees generated okay', {
  expect_true(all.equal(
    Preorder(ape::read.tree(text=("(0, (4, ((1, 5), (2, 3))));"))),
    as.phylo(10, 6, 0:5)))
  expect_equal(as.TreeNumber('10', 6, 0:5),
               as.TreeNumber(as.phylo.numeric(10, 6, 0:5)))
  Test <- function (i, nTip) {
    expect_equal(as.integer64(i),
                 as.integer64(as.TreeNumber(as.phylo(i, nTip, seq_len(nTip) - 1L))))
  }
  xx <- lapply(0:104, Test, 6)
  xx <- lapply(seq_len(NUnrooted(7)) - 1L, Test, 7)
  xx <- lapply(floor(runif(48) * NUnrooted(8)), Test, 8)
  xx <- lapply(floor(runif(48) * NUnrooted(10)), Test, 10)
  xx <- lapply(floor(runif(48) * NUnrooted(12)), Test, 12)

  nTip <- 14L
  treeNumber <- as.TreeNumber('123456789876', nTip, seq_len(nTip) - 1L)
  expect_equal(treeNumber, as.TreeNumber(as.phylo(treeNumber)))
  treeNumber <- as.TreeNumber('123456789876', nTip, seq_len(nTip) - 1L)
  expect_equal(treeNumber, as.TreeNumber(as.phylo(treeNumber)))
  
  expect_equal
})
