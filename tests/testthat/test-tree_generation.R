context('tree_generation.R')
test_that('Pectinate trees are generated', {
  expect_equal(ape::read.tree(text = '(t1, (t2, (t3, t4)));'),
               PectinateTree(4L))
  expect_equal(ape::read.tree(text = '(a, (b, (c, (d, e))));'),
               PectinateTree(letters[1:5]))
  expect_equal(ape::read.tree(text = '(a, (b, (c, (d, e))));'),
               PectinateTree(ape::read.tree(text = '(a, ((b, c), (d, e)));')))
  expect_equal(ape::read.tree(text = '(Cricocosmia, (Aysheaia, Siberion));'),
               PectinateTree(Lobo.phy[2:4]))
})

test_that('Balanced trees are generated correctly', {
  expect_equal(ape::read.tree(text = '(((t1, t2), (t3, t4)), ((t5,t6), (t7, (t8, t9))));'),
               BalancedTree(9L))
  expect_equal(integer(0), BalancedBit(seq_len(0)))
  expect_equal('Test', BalancedBit('Test'))
})
