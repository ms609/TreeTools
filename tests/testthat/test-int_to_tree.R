context('int_to_tree.cpp')

test_that('Trees generated okay', {
  expect_equal(as.phylo.numeric(10, 6, 0:5),
               Preorder(ape::read.tree(text=("(0, (4, ((1, 5), (2, 3))));"))))
  expect_equal(as.TreeNumber('10', 6, 0:5),
               as.TreeNumber(as.phylo.numeric(10, 6, 0:5)))
  Test <- function (i, nTip) {
    expect_equal(i, as.integer(as.TreeNumber(as.phylo(i, nTip, seq_len(nTip) - 1L))))
  }
  xx <- lapply(0:104, Test, 6)
  xx <- lapply(seq_len(NUnrooted(7)) - 1L, Test, 7)
  xx <- lapply(floor(runif(100) * NUnrooted(8)), Test, 8)
  xx <- lapply(floor(runif(100) * NUnrooted(10)), Test, 10)
  xx <- lapply(floor(runif(100) * NUnrooted(12)), Test, 12)

  nTip <- 14L
  treeNumber <- as.TreeNumber('123456789876', nTip, seq_len(nTip) - 1L)
  expect_equal(treeNumber, as.TreeNumber(as.phylo(treeNumber)))
})
