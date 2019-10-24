context('SplitFunctions.R')

test_that("SplitMatchProbability returns expected probabilities", {
  splitAB   <- c(rep(TRUE, 2), rep(FALSE, 7))
  splitABC  <- c(rep(TRUE, 3), rep(FALSE, 6))
  splitABI  <- c(rep(TRUE, 2), rep(FALSE, 6), TRUE)
  splitBCD  <- c(FALSE, rep(TRUE, 3), rep(FALSE, 5))
  splitAEF  <- c(TRUE, rep(FALSE, 3), rep(TRUE, 2), rep(FALSE, 3))
  splitABCD <- c(rep(TRUE, 4), rep(FALSE, 5))
  splitABCE <- c(rep(TRUE, 3), FALSE, TRUE, rep(FALSE, 4))
  splitCDEF <- c(rep(FALSE, 2), rep(TRUE, 4), rep(FALSE, 3))
  splitABEF <- c(rep(TRUE, 2), rep(FALSE, 2), rep(TRUE, 2), rep(FALSE, 3))
  splitABCDE<- c(rep(TRUE, 5), rep(FALSE, 4))

  splitAI <- c(TRUE, rep(FALSE, 7), TRUE)
  splitBC <- c(FALSE, TRUE, TRUE, rep(FALSE, 6))
  splitCD <- c(FALSE, FALSE, TRUE, TRUE, rep(FALSE, 5))

  expect_error(SplitMatchProbability(TRUE, splitAB))

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

  Test(1, splitAB, splitAI)
  Test(1, splitAB, splitBC)
  Test(1, splitBC, splitCD)
  Test(1, rev(splitAB), splitAI)
  Test(1L, splitABCD[-9], splitABEF[-9]) # Test even splits

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


  expect_equal(2L, MatchingSplitDistanceSplits(as.Splits(splitAB), as.Splits(splitAI)))
  expect_equal(2L, MatchingSplitDistanceSplits(splitAB, splitABCD))
  expect_equal(3L, MatchingSplitDistanceSplits(splitAB, splitABCDE))
  expect_equal(4L, MatchingSplitDistanceSplits(splitABC, splitAEF))
  expect_equal(MatchingSplitDistanceSplits(cbind(splitABC), cbind(splitAEF)),
               MatchingSplitDistanceSplits(cbind(splitAEF), cbind(splitABC)))
  expect_error(MatchingSplitDistanceSplits(cbind(splitAB), cbind(splitAB)[-9, ]))
  expect_error(NyeSplitSimilarity(cbind(splitAB), cbind(splitAB)[-9, ]))
})
