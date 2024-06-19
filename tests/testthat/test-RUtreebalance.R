test_that("Balance examples work as expected", {
  tree1 <- data.frame(Parent = c(1,1,1,1,2,3,4),
                      Identity = 1:7,
                      Population = c(1, rep(5, 6)))
  expect_equal(J1Index(tree1), 2/3)
  tree2 <- data.frame(Parent = c(1,1,1,1,2,3,4),
                      Identity = 1:7,
                      Population = c(rep(0, 4), rep(1, 3)))
  expect_equal(J1Index(tree2), 1/2)
  tree3 <- data.frame(Parent = c(1,1,1,1,2,3,4),
                      Identity = 1:7,
                      Population = c(0, rep(1, 3), rep(0, 3)))
  expect_equal(J1Index(tree3), 1)
  cat_tree <- data.frame(Parent = c(1, 1:14, 1:15, 15),
                         Identity = 1:31,
                         Population = c(rep(0, 15), rep(1, 16)))
  expect_equal(J1Index(cat_tree), 4736/9990)
  sym_tree <- data.frame(Parent = c(1, rep(1:15, each = 2)),
                         Identity = 1:31,
                         Population = c(rep(0, 15), rep(1, 16)))
  expect_equal(J1Index(sym_tree), 1)
  
  sym_tree2 <- data.frame(Parent = c(1, rep(1:15, each = 2)),
                         Identity = 1:31)
  expect_equal(J1Index(sym_tree), J1Index(sym_tree2))
  
})

test_that("Different tree inputs work", {
  phylo_tree2 <- read.tree(text='((A, B), ((C, D), (E, F)));')
  expect_equal(J1Index(phylo_tree2), 0.96936094)
  expected <- 0.79248125
  phylo_tree <- read.tree(text="((a:0.1)A:0.5,(b1:0.2,b2:0.1)B:0.2);")
  expect_equal(J1Index(phylo_tree), expected)
  phylo_tree <- read.tree(text="((1)5,(2,3)6);")
  expect_equal(J1Index(phylo_tree), expected)
  
  # data.frame omitting population sizes:
  edges_tree <- data.frame(Parent = c(4, 5, 4, 6, 6),
                           Identity = c(5, 1, 6, 2, 3))
  expect_equal(J1Index(edges_tree), expected)
  
  # data.frame omitting population sizes and including a row for the root:
  edges_tree_with_root <- data.frame(Parent = c(4, 4, 5, 4, 6, 6),
                                     Identity = c(4, 5, 1, 6, 2, 3))
  expect_equal(J1Index(edges_tree_with_root), expected)
  
  # data.frame including population sizes:
  edges_tree_with_pops <- data.frame(
    Parent = c(4, 5, 4, 6, 6),
    Identity = c(5, 1,6 ,2, 3),
    Population = c(0, 1, 0, 1, 1)
  )
  expect_equal(J1Index(edges_tree_with_pops), expected)
})

test_that(".GetAdjacency() works", {
  tree1 <- data.frame(Parent = c(1,1,1,1,2,3,4), Identity = 1:7,
                      Population = c(1, rep(5, 6)))
  expect_equal(
    .GetAdjacency(tree1),
    list(2:4, 5, 6, 7, NULL, NULL, NULL)
  )
})
