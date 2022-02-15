test_that('Subsplits', {
  splits <- as.Splits(PectinateTree(letters[1:9]))
  efgh <- Subsplit(splits, tips = letters[5:8], keepAll = TRUE, unique = FALSE)
  expect_equal(setNames(c(4, 4, 4, 3, 2, 1), 12:17),
               TipsInSplits(efgh))
  expect_equal(c('12' = TRUE, '13' = TRUE, '14' = TRUE, '15' = TRUE,
                 '16' = FALSE, '17' = TRUE), TrivialSplits(efgh))

  efghF <- Subsplit(splits, tips = letters[5:8], keepAll = FALSE)
  expect_equal(c('16' = 2), TipsInSplits(efghF))
  expect_equal(efghF, Subsplit(list(splits, splits), tips = letters[5:8],
                               keepAll = FALSE)[[1]])

  noSplit <- Subsplit(splits - splits, letters[5:8])
  expect_equal(attributes(noSplit)[3:5], attributes(efghF)[3:5])
  expect_equal(ignore_attr = TRUE, raw(1), noSplit[1])
  expect_equal(noSplit[1], Subsplit(splits - splits, 5:8)[1])


  splits <- as.Splits(PectinateTree(32 + 32 + 10))
  fourTips <- c('t32', 't33', 't64', 't65')
  sub <- Subsplit(splits, tips = fourTips)
  expect_equal(as.Splits(c('t32' = FALSE, 't33' = FALSE, 't64' = TRUE,
                           't65' = TRUE)),
               unname(sub), ignore_attr = TRUE)


})

test_that('Bitwise logic works', {
  splits <- as.Splits(BalancedTree(8))
  splits2 <- as.Splits(PectinateTree(8))
  A <- TRUE
  B <- FALSE
  .CSB <- function(a, b) .CompatibleSplit(packBits(as.logical(a)),
                                           packBits(as.logical(b)),
                                           length(a))
  expect_true(.CSB(rep(0, 8), rep(0, 8)))
  expect_true(.CSB(rep(0, 8), rep(1, 8)))
  expect_true(.CSB(rep(1, 8), rep(1, 8)))
  expect_true(.CSB(c(0, 0, 0, 0, 0, 0, 1, 1),
                   c(0, 0, 0, 0, 0, 0, 1, 1)))
  expect_true(.CSB(c(0, 0, 0, 0, 0, 0, 1, 1),
                   c(0, 0, 0, 0, 1, 1, 1, 1)))
  expect_true(.CSB(c(0, 0, 0, 0, 1, 1, 1, 1),
                   c(0, 0, 0, 0, 0, 0, 1, 1)))
  expect_true(.CSB(c(1, 1, 1, 1, 1, 0, 0, 0),
                   c(0, 0, 0, 1, 1, 1, 1, 1)))
  expect_false(.CSB(c(0, 0, 0, 0, 1, 1, 1, 1),
                    c(1, 0, 0, 0, 0, 0, 1, 1)))
  expect_false(.CSB(c(1, 0, 0, 0, 0, 0, 1, 1),
                    c(0, 0, 0, 0, 1, 1, 1, 1)))

  expectation <- matrix(TRUE, 5, 5,
                        dimnames = list(names(splits), names(splits2)))
  expectation["12", "12"] <- FALSE
  expectation["14", "14"] <- FALSE
  expect_equal(CompatibleSplits(splits, splits2), expectation)

  expect_true(.CompatibleSplit(as.raw(3), as.raw(7), nTip = 5))
  expect_false(.CompatibleSplit(as.raw(3), as.raw(6), nTip = 5))
})

test_that("SplitMatchProbability returns expected probabilities", {
  splitAB   <- as.Splits(c(rep(TRUE, 2), rep(FALSE, 7)))
  splitABC  <- as.Splits(c(rep(TRUE, 3), rep(FALSE, 6)))
  splitABI  <- as.Splits(c(rep(TRUE, 2), rep(FALSE, 6), TRUE))
  splitBCD  <- as.Splits(c(FALSE, rep(TRUE, 3), rep(FALSE, 5)))
  splitAEF  <- as.Splits(c(TRUE, rep(FALSE, 3), rep(TRUE, 2), rep(FALSE, 3)))
  splitABCD <- as.Splits(c(rep(TRUE, 4), rep(FALSE, 5)))
  splitABCE <- as.Splits(c(rep(TRUE, 3), FALSE, TRUE, rep(FALSE, 4)))
  splitCDEF <- as.Splits(c(rep(FALSE, 2), rep(TRUE, 4), rep(FALSE, 3)))
  splitABEF <- as.Splits(c(rep(TRUE, 2), rep(FALSE, 2), rep(TRUE, 2), rep(FALSE, 3)))
  splitABCDE<- as.Splits(c(rep(TRUE, 5), rep(FALSE, 4)))

  splitAI <- as.Splits(c(TRUE, rep(FALSE, 7), TRUE))
  splitHI   <- as.Splits(c(rep(FALSE, 7), rep(TRUE, 2)))
  splitBC <- as.Splits(c(FALSE, TRUE, TRUE, rep(FALSE, 6)))
  splitCD <- as.Splits(c(FALSE, FALSE, TRUE, TRUE, rep(FALSE, 5)))

  expect_error(SplitMatchProbability(as.Splits(TRUE), splitAB))

  expect_equal(LnSplitMatchProbability(splitABC, splitABCD),
               LnSplitMatchProbability(splitABCD,
                                       as.Splits(splitABC, rev(TipLabels(splitABC)))))

  # Possible matches to ABCD:....
  expect_true(SplitMatchProbability(splitABC, splitABCD) <
                SplitMatchProbability(splitAB, splitABCD))

  expect_true(SplitMatchProbability(splitABCD, splitABCD) <
                SplitMatchProbability(splitABC, splitABCD))

  expect_true(SplitMatchProbability(splitBCD, splitABCD) ==
                SplitMatchProbability(splitABC, splitABCD))

  expect_true(SplitMatchProbability(splitABC, splitABCD) <
                SplitMatchProbability(splitAEF, splitABCD))

  expect_true(SplitMatchProbability(splitABC, splitABCD) <
                SplitMatchProbability(splitABI, splitABCD))

  expect_true(SplitMatchProbability(splitABI, splitABCD) <
                SplitMatchProbability(splitAEF, splitABCD))

  expect_true(SplitMatchProbability(splitABC, splitABCD) <
                SplitMatchProbability(splitABCE, splitABCD))

  expect_true(SplitMatchProbability(splitCDEF, splitABCD) ==
                SplitMatchProbability(splitABEF, splitABCD))

  expect_true(SplitMatchProbability(splitABCE, splitABCD) <
                SplitMatchProbability(splitABEF, splitABCD))


  # Two splits of AB:...
  expect_true(SplitMatchProbability(splitAB, splitAB) <
                SplitMatchProbability(splitCD, splitAB))
  expect_true(SplitMatchProbability(splitCD, splitAB) <
                SplitMatchProbability(splitAI, splitAB))


  Test <- function(score, split1, split2) {
    expect_equal(ignore_attr = TRUE, score, SplitMatchProbability(split1, split2))
    expect_equal(ignore_attr = TRUE, score, SplitMatchProbability(split2, split1))

    expect_equal(ignore_attr = TRUE, score, SplitMatchProbability(split1, !split2))
    expect_equal(ignore_attr = TRUE, score, SplitMatchProbability(split2, !split1))

    expect_equal(ignore_attr = TRUE, score, SplitMatchProbability(!split1, !split2))
    expect_equal(ignore_attr = TRUE, score, SplitMatchProbability(!split2, !split1))

    expect_equal(ignore_attr = TRUE, score, SplitMatchProbability(!split1, split2))
    expect_equal(ignore_attr = TRUE, score, SplitMatchProbability(!split2, split1))

    score
  }

  Test(1L, splitAB, splitAI)
  Test(1L, splitAB, splitBC)
  Test(1L, splitBC, splitCD)
  Test(1L, splitHI, splitAI)
  Test(1L, as.Splits(c(rep(TRUE, 4), rep(FALSE, 4))), # Test even splits
       as.Splits(c(rep(TRUE, 2), rep(FALSE, 2), rep(TRUE, 2), rep(FALSE, 2))))

  Test(1/36, splitAB, splitAB)
  Test(1/36, splitBC, splitBC)

  Test(1/126, splitABCD, splitABCD)
  Test(1/126, splitABEF, splitABEF)
  Test(1/126, splitCDEF, splitCDEF)
  Test(1, splitABCD, splitABEF)

  Test(1/12, splitAB, splitABC)
  Test(1/12, splitBC, splitABC)
  Test(1/12, splitBC, splitBCD)

  Test(1, splitAEF, splitABCD)
  Test(66 / 126, splitABC, splitABEF)
  Test(66 / 126, splitBCD, splitCDEF)

  Test(4/84, splitABC, splitABCD)
  Test(4/84, splitBCD, splitABCD)
})


test_that('Tip labels are found', {
  pt4 <- PectinateTree(4L)
  bt4 <- BalancedTree(4L)
  t1..4 <- c('t1', 't2', 't3', 't4')
  expect_equal(t1..4, TipLabels(t1..4))
  expect_equal(t1..4, TipLabels(pt4))
  expect_equal(t1..4, TipLabels(as.Splits(pt4)))
  expect_equal(t1..4, TipLabels(list(pt4, bt4)))
  expect_equal(t1..4, TipLabels(structure(list(pt4, bt4),
                                          class='multiPhylo')))
  expect_equal(t1..4, TipLabels(structure(list(pt4, BalancedTree(letters[1:4])),
                                          class='multiPhylo'), single = TRUE))
  expect_equal(t1..4, TipLabels(structure(list(pt4, bt4,
                                               tip.label = t1..4),
                                          class='multiPhylo')))
  expect_equal(t1..4, TipLabels(as.TreeNumber(pt4)))

  expect_equal(t1..4, TipLabels(structure(list(), tip.label = t1..4)))
  expect_equal(t1..4, TipLabels(list(tip.label = t1..4)))
  expect_null(TipLabels(list()))
  expect_equal(TipLabels(pt4), TipLabels(list(as.Splits(pt4))))
  expect_equal(TipLabels(pt4, single = TRUE),
               TipLabels(list(as.Splits(pt4)), single = TRUE))
  expect_equal(TipLabels(pt4), TipLabels(list(as.Splits(pt4), as.Splits(bt4))))
  expect_equal(list(t1..4, c(t1..4, 't5')),
                    TipLabels(list(as.Splits(pt4),
                                   as.Splits(BalancedTree(5)))))

  expect_equal(t1..4, TipLabels(c(t1 = 1, t2 = 3, t3 = 3, t4 = 4)))
})

test_that("AllTipLabels()", {
  expect_equal(1:10, sort(as.integer(AllTipLabels(c(BalancedTree(1:5), PectinateTree(6:10))))))
  expect_equal(1:10, sort(as.integer(AllTipLabels(list(BalancedTree(1:6), PectinateTree(4:10))))))
})
