test_that("AddUnconstrained() works", {
  tips <- letters[1:9]
  constraint <- StringToPhyDat("0000?1111 000111111 0000??110", tips, FALSE)
  expect_equal(
    AddUnconstrained(constraint, letters[10:12], FALSE),
    PhyDatToMatrix(AddUnconstrained(constraint, letters[10:12], TRUE))
  )
  expect_equal(
    PhyDatToMatrix(AddUnconstrained(constraint, letters[8:12], TRUE)),
    AddUnconstrained(constraint, letters[10:12], FALSE)
  )
  
  noConstraint <- matrix(NA_character_, 6, 0,
                         dimnames = list(letters[1:6], NULL))
  expect_equal(
    AddUnconstrained(SingleTaxonTree("a"), letters[2:6], FALSE),
    noConstraint
  )
  
  expect_equal(AddUnconstrained(c(), letters[1:6], FALSE), noConstraint)
  expect_equal(AddUnconstrained(NULL, letters[1:6], FALSE), noConstraint)
  
  noNode <- structure(
    list(edge = structure(integer(0), dim = c(0L, 2L)),
         Nnode = 0L, tip.label = "a"), order = "preorder", class = "phylo")
  expect_equal(AddUnconstrained(noNode, letters[1:6], FALSE), noConstraint)
  
})

test_that("ImposeConstraint() works", {
  tips <- letters[1:9]
  tree <- as.phylo(1, 9, tips)
  
  expect_equal(ImposeConstraint(tree, c()), tree)
  expect_equal(ImposeConstraint(tree, KeepTip(tree, character(0))), tree)
  
  constraint <- StringToPhyDat("0000?1111 000111111 0000??110", tips, FALSE)
  expect_true(all.equal(
    ImposeConstraint(tree, constraint),
    read.tree(text = "((a, (b, c)), (d, (e, (f, (i, (g, h))))));")))

  expect_equal(ImposeConstraint(tree, constraint),
               ImposeConstraint(tree, PhyDatToMatrix(constraint)))
  
  expect_equal(
    ImposeConstraint(tree, setNames(c(rep(0, 3), rep(1, 6)), tips)),
    ImposeConstraint(tree, StringToPhyDat("000111111", tips, FALSE))
  )

  constraint <- StringToPhyDat("00001111 00011111 0000?110", tips[-5], FALSE)
  expect_true(all.equal(
    ImposeConstraint(tree, constraint),
    read.tree(text = "((a, (b, c)), (d, (e, (f, (i, (g, h))))));")))

  constraint <- ape::read.tree(text = "(a, b, c, (d, (f, (g, h), i)));")
  expect_true(all.equal(
    ImposeConstraint(tree, constraint),
    read.tree(text = "((a, (b, c)), (d, (e, (f, (i, (g, h))))));")))
  
  constraint <- MatrixToPhyDat(matrix(
    c(0, 0, "?", "?", 1, 1,
      1, 1,   1, "?", 0, 0), ncol = 2,
    dimnames = list(letters[1:6], NULL)))

  expect_warning({out <- ImposeConstraint(NJTree(constraint), constraint)})
  expect_true(inherits(out, "phylo"))

  # Need to collapse splits intelligently to avoid error.  From Joe Moysiuk.
  constraint <- MatrixToPhyDat(
    structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L,
                1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 1L,
                1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L,
                1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L,
                1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
              .Dim = c(38L, 15L),
              .Dimnames = list(c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8",
                                 "t9", "t10", "t11", "t12", "t13", "t14", "t15",
                                 "t16", "t17", "t18", "t19", "t20", "t21",
                                 "t22", "t23", "t24", "t25", "t26", "t27",
                                 "t28", "t29", "t30", "t31", "t32", "t33",
                                 "t34", "t35", "t36", "t37", "t38"), NULL)))
  ImposeConstraint(NJTree(constraint), constraint)
})
