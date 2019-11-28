context('int_to_tree.cpp')

test_that('Trees generated okay', {
  expect_equal(as.phylo.numeric(10, 6, 0:5),
               Preorder(ape::read.tree(text=("(0, (4, ((1, 5), (2, 3))));"))))
})
