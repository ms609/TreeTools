context('tree_generation.R')
test_that('Pectinate trees are generated', {
  expect_equal(ape::read.tree(text = '(t1, (t2, (t3, t4)));'),
               PectinateTree(4L))
  expect_equal(ape::read.tree(text = '(a, (b, (c, (d, e))));'),
               PectinateTree(letters[1:5]))
  expect_equal(ape::read.tree(text = '(Cricocosmia, (Aysheaia, Siberion));'),
               PectinateTree(Lobo.phy[2:4]))
})