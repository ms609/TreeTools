context('SplitFunctions.R')

test_that('Subsplits', {
  splits <- as.Splits(PectinateTree(letters[1:9]))
  efgh <- Subsplit(splits, tips = letters[5:8], keepAll = TRUE)
  expect_equal(c(4, 4, 4, 3, 2, 1), as.integer(TipsInSplits(efgh)))
  expect_equal(c(n12 = TRUE, n13 = TRUE, n14 = TRUE, n15 = TRUE, n16 = FALSE,
                 n17 = TRUE), TrivialSplits(efgh))

  efghF <- Subsplit(splits, tips = letters[5:8], keepAll = FALSE)
  expect_equal(c(n16 = 2), TipsInSplits(efghF))
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
