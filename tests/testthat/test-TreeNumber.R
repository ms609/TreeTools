context("TreeNumber.R")

test_that("as.phylo.numeric()", {
  expect_equal(as.phylo(0:2, 6, letters[1:6])[[1]],
               PectinateTree(letters[c(1, 3:6, 2)]))
  expect_equal(as.phylo(0, tipLabels = letters[1:6]),
               PectinateTree(letters[c(1, 3:6, 2)]))
  expect_error(as.phylo(0))
})

test_that("as.TreeNumber()", {
  expect_equal(c("Phylogenetic tree number 0 of 105 ",
                 " 6 tips: t1 t2 t3 t4 t5 t6"),
               capture.output(print(as.TreeNumber(as.phylo(105, 6)))))
  expect_equal(1:3, unlist(as.TreeNumber(as.phylo(1:3, 6))))
  tn <- as.TreeNumber(as.phylo(1, 6))
  attr(tn, 'tip.label') <- NULL
  expect_equal(as.phylo(1, 6), as.phylo(tn))
})
