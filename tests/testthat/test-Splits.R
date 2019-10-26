context("Splits.R")

test_that("as.Split", {
  A <- FALSE
  B <- TRUE
  expect_equal(c("1 bipartition split dividing 4 tips.",
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
  expect_equal(letters[1:5], colnames(as.logical(logicalSplits)))
})

test_that("Split operations", {
  split1 <- as.Splits(c(rep(TRUE, 30), rep(FALSE, 70), rep(TRUE, 20)))
  notSplit1 <- as.Splits(c(rep(FALSE, 30), rep(TRUE, 70), rep(FALSE, 20)))
  expect_equal(!split1,notSplit1)
  expect_equal(split1,!notSplit1)
  
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
  expect_error(c(splits1, as.Splits(tree3)))
  expect_error(c(splits1, as.Splits(tree4)))
  expect_error(c(splits1, as.Splits(tree5)))
  expect_equal(c(28L, 24L, 28L, 24L),
               as.integer(c(splits1,
                            as.Splits(RenumberTips(tree3, letters[1:5])))[, 1]))
  expect_equal(2L, length(unique(c(splits1, as.Splits(tree2)))))
  expect_equal(2L, length(unique(c(splits1,
                                   as.Splits(c(3, 7), tipLabels=letters[1:5])))))


  # TODO: Fully test splits with large (>32 tip) trees
})
