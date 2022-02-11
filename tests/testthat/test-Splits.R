expect_splits_equal <- function(s1, s2, ...) {
  expect_equal(as.character(PolarizeSplits(s1)),
               as.character(PolarizeSplits(s2)), ...)
}

test_that("as.Splits()", {
  A <- FALSE
  B <- TRUE
  expect_equal(strsplit(capture_output(summary(as.Splits(c(A, A, B, B)))),
                        '\n')[[1]],
               c("1 bipartition split dividing 4 tips, t1 .. t4",
                 "   1234",
                 "   ..**", "",
                 " Tip 1: t1\t Tip 2: t2\t Tip 3: t3\t Tip 4: t4\t"))
  logical80 <- c(rep(TRUE, 40), rep(FALSE, 16), rep(TRUE, 24))
  expect_equal(strsplit(capture_output(print(
    as.Splits(logical80), detail = TRUE)), '\n')[[1]][3],
    paste0(c('   ', ifelse(logical80, '*', '.')), collapse = '')
  )
  expect_equal(as.logical(as.logical(as.Splits(logical80))), logical80)
  expect_equal(as.logical(as.Splits(c(A, A, B, B))),
               t(matrix(c(A, A, B, B), dimnames = list(paste0('t', 1:4)))))
  tree1 <- BalancedTree(letters[1:5])
  splits1 <- as.Splits(tree1)
  expect_equal(as.character(splits1),
               c("7" = 'a b c | d e', "8" = "a b | c d e"))
  
  logicalSplits <- as.Splits(matrix(c(B, B, B, A, A,  B, B, A, A, A),
                                    nrow = 2, byrow = TRUE),
                             tipLabels = letters[1:5])
  rownames(logicalSplits) <- rownames(splits1)
  expect_equal(splits1, logicalSplits)
  expect_equal(splits1, as.Splits(splits1))
  
  pec6 <- as.Splits(PectinateTree(letters[1:6]))
  expect_equal(as.character(as.Splits(pec6)),
               c("9" = "c d e f | a b",
                 "10" = "d e f | a b c",
                 "11" = "e f | a b c d"))
  
  splitsC <- as.Splits(ape::read.tree(text = "(((a, d), e), (b, (f, c)));"))
  splitsD <- as.Splits(ape::read.tree(text = "((a, b, c), (d, (e, f)));"))
  splitsU <- as.Splits(ape::read.tree(text = "(a, b, c, d, e, f);"))
  oneSplit <- as.Splits(ape::read.tree(text = "((a, b, c), (d, e, f));"))
  expect_equal(attr(splitsC, 'tip.label'),
               attr(as.Splits(splitsU, splitsC), 'tip.label'))
  expect_equal(attr(splitsC, 'tip.label'),
               attr(as.Splits(oneSplit, splitsC), 'tip.label'))
  expect_equal(as.Splits(splitsD, splitsC),
               as.Splits(list(splitsC, splitsD, splitsU, oneSplit))[[2]])

  expect_equal(letters[1:5], colnames(as.logical(logicalSplits)))

  polytomy <- ape::read.tree(text = '(a, b, c, d, e);')
  expect_equal("0 bipartition splits dividing 5 tips, a .. e",
               capture_output(print(as.Splits(polytomy))))

  notPreorder <- structure(list(
    edge = structure(c(6L, 9L, 8L, 7L, 7L, 8L, 9L, 6L,
                       9L, 8L, 7L, 2L, 3L, 5L, 4L, 1L), .Dim = c(8L, 2L)),
    Nnode = 4L, tip.label = 1:5), class = "phylo", order = "cladewise")
  expect_equal(c('7' = packBits(c(A, B, B, A, A, rep(FALSE, 3))),
                 '8' = packBits(c(A, B, B, A, B, rep(FALSE, 3)))),
               as.Splits(notPreorder)[, 1])
})

test_that('as.Splits.phylo()', {
  rootedStar <- structure(list(edge = structure(c(7L, 8L, 8L, 8L, 8L, 8L, 7L,
                                                  8L, 1L, 2L, 3L, 4L, 5L, 6L),
                                                .Dim = c(7L, 2L)), Nnode = 2L,
                               tip.label = letters[1:6]), class = "phylo")

  rootedStar2 <- structure(list(edge = structure(c(8L, 8L, 8L, 8L, 8L, 7L, 7L,
                                                   1L, 2L, 3L, 4L, 5L, 6L, 8L),
                                                .Dim = c(7L, 2L)), Nnode = 2L,
                               tip.label = letters[1:6]), class = "phylo")
  nasty <- structure(list(edge = structure( # Danger: Do not plot!
    c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
      5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
    .Dim = c(12, 2)), Nnode = 5L, tip.label = letters[1:8]), class = 'phylo')


  expect_equal(c(0L, 1L), dim(as.Splits(rootedStar)))
  expect_equal(c(0L, 1L), dim(as.Splits(rootedStar2)))
  expected <- as.raw(as.Splits(Preorder(nasty)))
  actual <- as.raw(as.Splits(nasty))
  expect_equal(length(expected), length(actual))
  expect_equal(length(actual), length(unique(actual)))
  expect_true(all(actual %in% expected))
  expect_true(all(expected %in% actual))

  expect_equal(c(2L, 1L), dim(as.Splits(PectinateTree(5L))))
  expect_equal(c(2L, 1L), dim(as.Splits(UnrootTree(PectinateTree(5L)))))

  twoCherriesInGrass <- PectinateTree(3) + PectinateTree(3)
  expect_equal(c(2L, 1L), dim(as.Splits(twoCherriesInGrass)))

  unresolvedRoot <- CollapseNode(twoCherriesInGrass, 8)
  expect_equal(c(1L, 1L), dim(as.Splits(unresolvedRoot)))
  expect_equal(c(3L, 1L), dim(as.Splits(BalancedTree(6))))


  expect_equal(c(61L, 8L), dim(as.Splits(PectinateTree(64L))))
  expect_equal(c(62L, 9L), dim(as.Splits(PectinateTree(65L))))
  expect_equal(c(125L, 16L), dim(as.Splits(PectinateTree(128L))))
  expect_equal(as.Splits(Postorder(PectinateTree(64L))),
               as.Splits(PectinateTree(64L)))
  expect_true(all(
    as.Splits(ape::reorder.phylo(BalancedTree(66L), 'postorder')) %in% 
               as.Splits(BalancedTree(66L))))
  expect_false(all(
    as.Splits(ape::reorder.phylo(PectinateTree(66L), 'postorder')) %in% 
               as.Splits(BalancedTree(66L))))

  # Should be viable on 32-bit systems
  expect_equal(c(7997L, 1000L), dim(as.Splits(PectinateTree(8000L))))

  bigTree <- PectinateTree(16384L) # Max size with int16_t on 64-bit systems
  bigMemory <- try(Postorder(bigTree$edge), silent = TRUE)
  if (inherits(bigMemory, 'try-error')) {
    expect_error(dim(as.Splits(bigTree))) # Likely on 32-bit systems
  } else {
    rm(bigMemory)
    gc()
    expect_equal(c(16381L, 2048L), dim(as.Splits(bigTree)))
  }

})

test_that('as.Splits.multiPhylo()', {
  randomTreeIds <- c(30899669, 9149275, 12823175, 19740197, 31296318,
                     6949843, 30957991, 32552966, 22770711, 21678908)
  randomTrees <- as.phylo(randomTreeIds, 11L, seq_len(11L))

  expect_equal(as.Splits(randomTrees[[1]]),
               as.Splits(randomTrees)[[1]])
})

test_that('as.Splits.Splits()', {
  # n tip > 8L
  splitsA <- as.Splits(ape::read.tree(text="((((a, b, c, c2), g), h), (d, (e, f)));"))
  splitsB <- as.Splits(ape::read.tree(text="(((((a, b), (c, c2)), h), g), (d, e, f));"))
  expectedSplits <- structure(matrix(as.raw(c(0x3f, 0x2f, 0x0f, 0x03, 0x0c, 0,
                                              0, 0, 0, 0)), ncol=2), nTip = 9,
                              tip.labels = TipLabels(splitsA), class='Splits')
  actualSplits <- as.Splits(splitsB, splitsA)
  expect_true(all(expectedSplits %in% actualSplits))
  expect_true(all(actualSplits %in% expectedSplits))
  expect_equal(NSplits(expectedSplits), NSplits(actualSplits))

  splitsABare <- splitsA
  attr(splitsABare, 'tip.label') <- NULL
  expect_error(as.Splits(splitsABare, 1:2))
  expect_equal(splitsA, as.Splits(splitsABare, splitsA))
})

test_that('as.Splits.matrix()', {
  expect_error(as.Splits(matrix(1, 3, 3)))
  expect_error(as.Splits(matrix(1, 2, 2,
                                dimnames = list(c('edge', 'Nnode'), NULL))))
  trees <- list(BalancedTree(8), PectinateTree(8),
                CollapseNode(BalancedTree(8), 10:13))
  arr <- sapply(1:3, function(i) trees[[i]])

  expect_error(as.Splits(arr[-2, ]))
  expect_equal(as.Splits(trees), as.Splits(arr))

  for (i in seq_along(trees)) {
    trees[[i]]$tip.label <- 1:8
  }
  expect_equal(as.Splits(trees), as.Splits(arr[-3, ]))

  exp <- as.Splits(BalancedTree(8))
  attr(exp, 'tip.label') <- NULL
  expect_equal("Tips not labelled.", capture.output(summary(exp))[9])
  expect_equal(exp, as.Splits(BalancedTree(8)$edge))
})

test_that('as.Splits.edge()', {
  test <- unname(as.Splits(BalancedTree(4), asSplits = FALSE))
  expect_true(identical(test, matrix(as.raw(0x0c))) ||
              identical(test, matrix(as.raw(0x03))))
})

test_that('as.Splits.logical()', {
  FFTT <- c(FALSE, FALSE, TRUE, TRUE)
  a..d <- letters[1:4]
  expect_splits_equal(unname(as.Splits(BalancedTree(a..d))),
                      as.Splits(!FFTT, a..d))
  expect_splits_equal(as.Splits(FFTT, a..d), as.Splits(t(matrix(FFTT)), a..d))
  expect_splits_equal(as.Splits(FFTT), as.Splits(t(matrix(FFTT))))
})

test_that("&.Splits()", {
  splits <- structure(as.raw(c(0x07, 0x03, 0x18, 0xe0, 0x60, 0x80)),
                      .Dim = c(3L, 2L), nTip = 9L, class = "Splits")
  mask <- as.raw(c(0x0f, 0x00))
  expect_equal(structure(as.raw(c(0x07, 0x00, 0x03, 0x00, 0x08, 0x00)),
                         .Dim = c(2L, 3L)),
               t(splits) & mask)
})

test_that('empty as.X.Splits()', {
  someSplit <- as.Splits(BalancedTree(6))
  noSplit <- someSplit[[logical(0)]]
  expect_equal(matrix(logical(0), 0, 6,
                      dimnames = list(NULL, paste0('t', 1:6))),
               as.logical(noSplit))
  expect_equal(character(0), as.character(noSplit))
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

test_that('match.Splits()', {
  tree1 <- PectinateTree(1:8)
  tree2 <- PectinateTree(8:1)
  col2 <- as.Splits(CollapseNode(tree2, 13))

  expect_equal(5:1, match(as.Splits(tree1), as.Splits(tree2, tree1)))
  expect_equal(c(4, 3, NA, 2, 1), match(as.Splits(tree1, tree2), col2))
  expect_equal(c(4, 3, 999, 2, 1),
               match(as.Splits(tree1, tree2), col2, nomatch = 999))
  expect_equal(c(5, 4, 2, 1), match(col2, as.Splits(tree1, tree2)))
})

test_that("%in%.Splits()", {
  splits5 <- as.Splits(PectinateTree(5))
  x <- splits5
  table <- splits5[[-1]]
  expect_null(dim(duplicated(c(splits5, splits5[[-1]]))))
  expect_equal(names(splits5 %in% splits5[[-1]]), names(splits5))
  
  splits9 <- as.Splits(PectinateTree(9))
  expect_equal(names(splits9 %in% splits9[[-4]]), names(splits9))
})

test_that("print.Splits()", {
  sp4 <- as.Splits(BalancedTree(4))
  num <- names(sp4)
  expect_equal(
    capture.output(print(PolarizeSplits(sp4, 1), details = TRUE)),
    c( "1 bipartition split dividing 4 tips, t1 .. t4", "    1234",
       paste0(" ", num, "  **..")))
  expect_equal(capture.output(print(PolarizeSplits(sp4, 4), details = TRUE)),
               c( "1 bipartition split dividing 4 tips, t1 .. t4", "    1234",
                  paste0(" ", num, "  ..**")))
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
  if (isTRUE(all.equal(rownames(namedSplits), 
                       as.character(c(10:12, 14:15))))) {
    expect_true(TRUE)
  } else {
    expect_equal(rownames(namedSplits), as.character(c(11:15)))
  }
  
  expect_equal(split1[], `[[.Splits`(split1))
  expect_error(c(split1, notSplit1)[[2:1, 1]])
})

test_that("Split subtraction", {
  splits <- as.Splits(BalancedTree(8))
  expect_equal(splits[[c(1, 2, 4, 5)]], splits - splits[[3]])
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
  expect_equal(c(as.raw(splits1), as.raw(as.Splits(PectinateTree(5)))),
               unname(c(splits1,
                        as.Splits(RenumberTips(tree3, letters[1:5])))[, 1]))
  expect_equal(2L, length(unique(c(splits1, as.Splits(tree2)))))

  #TODO: Fully test splits with large (> 8 tip) trees
})

test_that("as.phylo.Splits()", {
  Test <- function(tr) {
    expect_true(all.equal(as.phylo(as.Splits(tr)), unroot(tr)))
  }
  Test(read.tree(text = "((a, b, c), (d, e, f, g));"))
  Test(read.tree(text = "((a, b, c), (d, e, (f, g)));"))
  Test(BalancedTree(64))
  Test(BalancedTree(65)) # multiple bins
  Test(PectinateTree(65*6)) # big!

  A <- FALSE
  B <- TRUE
  splits <- as.Splits(matrix(c(A, B, A, A, A, A, A, B,
                               A, B, A, A, A, A, B, B,
                               A, B, A, A, A, B, B, B,
                               A, B, A, A, B, B, B, B), 4, byrow = TRUE),
                      tipLabels = letters[1:8])
  expect_true(all.equal(as.phylo(splits),
                        read.tree(text = '(a, c, d, (e, (f, (g, (h, b)))));')))
})

test_that("PolarizeSplits()", {
  bal6 <- as.Splits(BalancedTree(6))
  expect_error(PolarizeSplits(bal6, 'ERROR'))
  expect_error(PolarizeSplits(bal6, 0))
  expect_error(PolarizeSplits(bal6, 7))
})
