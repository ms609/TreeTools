context("Splits.R")

test_that("as.Split", {
  A <- FALSE
  B <- TRUE
  expect_equal(c("1 bipartition split dividing 4 tips, t1 .. t4",
                 "   1234",
                 "   ..**", "",
                 " Tip 1: t1\t Tip 2: t2\t Tip 3: t3\t Tip 4: t4\t"),
               strsplit(capture_output(summary(as.Splits(c(A, A, B, B)))), '\n')[[1]])
  logical80 <- c(rep(TRUE, 40), rep(FALSE, 16), rep(TRUE, 24))
  expect_equal(
    paste0(c('   ', ifelse(logical80, '*', '.')), collapse=''),
    strsplit(capture_output(print(as.Splits(logical80), detail = TRUE)), '\n')[[1]][3]
  )
  expect_equal(as.logical(as.logical(as.Splits(logical80))), logical80)
  expect_equal(t(matrix(c(A, A, B, B), dimnames=list(paste0('t', 1:4)))),
               as.logical(as.Splits(c(A, A, B, B))))
  tree1 <- BalancedTree(letters[1:5])
  splits1 <- as.Splits(tree1)
  expect_equal(c(n8 = 'c d e | a b', n9 = 'd e | a b c'), as.character(splits1))
  logicalSplits <- as.Splits(matrix(c(A, A, B, B, B,  A, A, A, B, B),
                                    nrow=2, byrow = TRUE),
                             tipLabels = letters[1:5])
  rownames(logicalSplits) <- rownames(splits1)
  expect_equal(splits1, logicalSplits)
  expect_equal(splits1, as.Splits(splits1))

  splitsC <- as.Splits(ape::read.tree(text="(((a, d), e), (b, (f, c)));"))
  splitsD <- as.Splits(ape::read.tree(text="((a, b, c), (d, (e, f)));"))
  splitsU <- as.Splits(ape::read.tree(text="(a, b, c, d, e, f);"))
  oneSplit <- as.Splits(ape::read.tree(text="((a, b, c), (d, e, f));"))
  expect_equal(attr(splitsC, 'tip.label'),
               attr(as.Splits(splitsU, splitsC), 'tip.label'))
  expect_equal(attr(splitsC, 'tip.label'),
               attr(as.Splits(oneSplit, splitsC), 'tip.label'))

  expect_equal(letters[1:5], colnames(as.logical(logicalSplits)))

  polytomy <- ape::read.tree(text='(a, b, c, d, e);')
  expect_equal("0 bipartition splits dividing 5 tips, a .. e",
               capture_output(print(as.Splits(polytomy))))

  notPreOrder <- structure(list(edge = structure(c(6L, 9L, 8L, 7L, 7L, 8L, 9L,
                                            6L, 9L, 8L, 7L, 2L, 3L, 5L, 4L, 1L),
                                          .Dim = c(8L, 2L)), Nnode = 4L,
                         tip.label = 1:5), class = "phylo", order = "cladewise")
  expect_equal(c(n8 = 22, n9 = 6), as.Splits(notPreOrder)[, 1])

  expect_equal(c(61L, 2L), dim(as.Splits(PectinateTree(64L))))
  expect_equal(c(125L, 4L), dim(as.Splits(PectinateTree(128L))))
})


test_that('Renaming splits', {
  tree1 <- PectinateTree(1:8)
  tree2 <- BalancedTree(8:1)
  splits1 <- as.Splits(tree1)
  splits2 <- as.Splits(tree2)

  expect_error(as.Splits(splits1, tipLabel = 1:4))
  expect_error(as.Splits(splits1, tipLabel = LETTERS[1:8]))

  expect_equal(as.Splits(tree1, tipLabel = as.character(8:1)),
               as.Splits(tree1, tipLabel = tree2))

  expect_equal(as.Splits(tree1, tipLabel = tree2),
               as.Splits(as.Splits(tree1), tipLabel = tree2))

  as.Splits(tree1, tipLabel = tree2)[]
  as.Splits(as.Splits(tree1), tipLabel = tree2)[]

})

test_that('match.Splits', {
  tree1 <- BalancedTree(1:8)
  tree2 <- BalancedTree(8:1)
  col2 <- as.Splits(CollapseNode(tree2, 13))

  expect_equal(5:1, match.Splits(as.Splits(tree1), as.Splits(tree2, tree1)))
  expect_equal(c(4, 3, NA, 2, 1), match.Splits(as.Splits(tree1, tree2), col2))
  expect_equal(c(5, 4, 2, 1), match.Splits(col2, as.Splits(tree1, tree2)))
})

test_that("Split operations", {
  split1 <- as.Splits(c(rep(TRUE, 30), rep(FALSE, 70), rep(TRUE, 20)))
  notSplit1 <- as.Splits(c(rep(FALSE, 30), rep(TRUE, 70), rep(FALSE, 20)))
  expect_equal(!split1, notSplit1)
  expect_equal(split1, !notSplit1)

  expect_equal(notSplit1[1, ], (split1 + notSplit1)[2, ])
  expect_equal(split1 + notSplit1, !(notSplit1 + split1))

  expect_equal(notSplit1, (split1 + notSplit1)[[2]])
  expect_equal(notSplit1 + split1 + split1, (split1 + notSplit1)[[c(2, 1, 1)]])

  split2 <- as.Splits(c(rep(TRUE, 4), FALSE, FALSE))
  notSplit2 <- as.Splits(c(rep(FALSE, 4), TRUE, TRUE))

  expect_equal(split2, !notSplit2)
  expect_equal(!split2, notSplit2)
  expect_equal(notSplit2[1, ], (split2 + notSplit2)[2, ])
  expect_equal(split2 + notSplit2, !(notSplit2 + split2))

  expect_error(split1 + split2)

  namedSplits <- as.Splits(BalancedTree(8))
  expect_equal(c('n12', 'n14'), rownames(namedSplits[[c(2, 4)]]))

})

test_that("Split combination", {
  tree1 <- BalancedTree(letters[1:5])
  splits1 <- as.Splits(tree1)
  tree2 <- PectinateTree(letters[1:5])
  splits12 <- c(splits1, as.Splits(tree2))
  tree3 <- PectinateTree(letters[c(1:3, 5, 4)])
  tree4 <- PectinateTree(letters[c(1:4, 6)])
  tree5 <- PectinateTree(letters[1:6])

  expect_equal(4L, length(splits12))
  expect_equal(c(FALSE, FALSE, TRUE, TRUE), as.logical(duplicated(splits12)))
  expect_equal(2L, length(unique(c(splits1, as.Splits(tree3)))))
  expect_error(c(splits1, as.Splits(tree4)))
  expect_error(c(splits1, as.Splits(tree5)))
  expect_equal(c(28L, 24L, 28L, 24L),
               as.integer(c(splits1,
                            as.Splits(RenumberTips(tree3, letters[1:5])))[, 1]))
  expect_equal(2L, length(unique(c(splits1, as.Splits(tree2)))))
  expect_equal(2L, length(unique(c(splits1,
                                   as.Splits(c(3, 7), tipLabels=letters[1:5])))))

  expect_equal(c(n8 = 3, n9 = 2), TipsInSplits(splits1))

  # TODO: Fully test splits with large (>32 tip) trees
})
