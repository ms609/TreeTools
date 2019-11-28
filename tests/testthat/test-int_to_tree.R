context('int_to_tree.cpp')

test_that('Trees generated okay', {
  expect_equal(as.phylo.numeric(10, 6, 0:5),
               Preorder(ape::read.tree(text=("(0, (4, ((1, 5), (2, 3))));"))))
  Test <- function (i, nTip) {
    expect_equal(i, as.integer(as.numeric.phylo(as.phylo.numeric(i, nTip, seq_len(nTip) - 1L))))
  }
  lapply(0:104, Test, 6)
  lapply(seq_len(NUnrooted(7)) - 1L, Test, 7)

})
