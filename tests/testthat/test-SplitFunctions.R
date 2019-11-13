context('SplitFunctions.R')

test_that('Subsplits', {
  splits <- as.Splits(PectinateTree(letters[1:9]))
  efgh <- Subsplit(splits, tips = letters[5:8], keepAll = TRUE, unique = FALSE)
  expect_equal(c(4, 4, 4, 3, 2, 1), as.integer(TipsInSplits(efgh)))
  expect_equal(c(n12 = TRUE, n13 = TRUE, n14 = TRUE, n15 = TRUE, n16 = FALSE,
                 n17 = TRUE), TrivialSplits(efgh))

  efghF <- Subsplit(splits, tips = letters[5:8], keepAll = FALSE)
  expect_equal(c(n16 = 2), TipsInSplits(efghF))

  splits <- as.Splits(PectinateTree(32+32+10))
  sub <- Subsplit(splits, tips = c('t32', 't33', 't64', 't65'))
})

test_that('Bitwise logic works', {
  splits <- as.Splits(BalancedTree(8))
  splits2 <- as.Splits(PectinateTree(8))
  A <- TRUE
  B <- FALSE
  dn <- paste0('n', 11:15)
  .AsDecimal <- function (bool) sum(2^(8:0) * bool)
  .CSB <- function (a, b) .CompatibleSplit(.AsDecimal(a), .AsDecimal(b),
                                           length(a))
  expect_true(.CSB(rep(0, 9), rep(0, 9)))
  expect_true(.CSB(rep(0, 9), rep(1, 9)))
  expect_true(.CSB(rep(1, 9), rep(1, 9)))
  expect_true(.CSB(c(0, 0, 0, 0, 0, 0, 0, 1, 1),
                   c(0, 0, 0, 0, 0, 0, 0, 1, 1)))
  expect_true(.CSB(c(0, 0, 0, 0, 0, 0, 0, 1, 1),
                   c(0, 0, 0, 0, 0, 1, 1, 1, 1)))
  expect_true(.CSB(c(0, 0, 0, 0, 0, 1, 1, 1, 1),
                   c(0, 0, 0, 0, 0, 0, 0, 1, 1)))
  expect_true(.CSB(c(1, 1, 1, 1, 1, 0, 0, 0, 0),
                   c(0, 0, 0, 1, 1, 1, 1, 1, 1)))
  expect_false(.CSB(c(0, 0, 0, 0, 0, 1, 1, 1, 1),
                    c(1, 0, 0, 0, 0, 0, 0, 1, 1)))
  expect_false(.CSB(c(1, 0, 0, 0, 0, 0, 0, 1, 1),
                    c(0, 0, 0, 0, 0, 1, 1, 1, 1)))

  expect_equal(
    matrix(c(A, A, A, A, A,
             A, B, A, A, A,
             A, A, A, A, A,
             A, A, A, B, A,
             A, A, A, A, A), byrow = TRUE, 5, 5, dimnames = list(dn, dn)),
    CompatibleSplits(splits, splits2))
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


  Test <- function (score, split1, split2) {
    expect_equivalent(score, SplitMatchProbability(split1, split2))
    expect_equivalent(score, SplitMatchProbability(split2, split1))

    expect_equivalent(score, SplitMatchProbability(split1, !split2))
    expect_equivalent(score, SplitMatchProbability(split2, !split1))

    expect_equivalent(score, SplitMatchProbability(!split1, !split2))
    expect_equivalent(score, SplitMatchProbability(!split2, !split1))

    expect_equivalent(score, SplitMatchProbability(!split1, split2))
    expect_equivalent(score, SplitMatchProbability(!split2, split1))

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
  expect_equal(c('t1', 't2', 't3', 't4'), TipLabels(c('t1', 't2', 't3', 't4')))
  expect_equal(c('t1', 't2', 't3', 't4'), TipLabels(pt4))
  expect_equal(c('t1', 't2', 't3', 't4'), TipLabels(as.Splits(pt4)))
  expect_equal(c('t1', 't2', 't3', 't4'), TipLabels(list(pt4, bt4)))
  expect_equal(c('t1', 't2', 't3', 't4'), TipLabels(structure(list(pt4, bt4),
                                                    class='multiPhylo')))
  expect_equal(TipLabels(pt4), TipLabels(list(as.Splits(pt4))))
})
