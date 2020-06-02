context("combinatorics.R")

test_that("Factorials are calculated correctly", {
  expect_equal(c(1L, 1L, 1L, 2L, 3L, 2L * 4L, 3L * 5L,
                 2L * 4L * 6L, 3L * 5L * 7L), DoubleFactorial(-1:7))
  expect_equal(doubleFactorials[1:20], DoubleFactorial(1:20))
  expect_equal(LnDoubleFactorial(-1:10), log(DoubleFactorial(-1:10)))
  expect_equal(logDoubleFactorials[1:20], log(DoubleFactorial(1:20)))
  expect_equal(log2DoubleFactorials[1:20], log2(DoubleFactorial(1:20)))
  expect_equal(LnDoubleFactorial(50001) - log(50001), LnDoubleFactorial.int(49999L))
  expect_equal(Log2DoubleFactorial(50000 - 2 + -1:1) + log2(50000 + -1:1),
               Log2DoubleFactorial(50000 + -1:1))
  expect_equal(LnDoubleFactorial.int(-1L), 0L)
  expect_equal(LnRooted.int(-1L), 0L)

  expect_error(DoubleFactorial(301))
})

test_that("N consistent with splits calculated correctly", {
  expect_equal(NUnrooted(8), NUnrootedSplits(8))
  expect_equal(NRooted(3) * NRooted(5), NUnrootedSplits(3, 5))
  expect_equal(doubleFactorials[12 + 12 - 5] /
    doubleFactorials[17] * doubleFactorials[4 + 4 - 3] ^ 3,
         NUnrootedSplits(c(4, 4, 4)))

  OldNUS <- function (...) {
    splits <- c(...)
    splits <- splits[splits > 0]
    totalTips <- sum(splits)
    round(DoubleFactorial(totalTips + totalTips - 5L) /
            DoubleFactorial(2L * (totalTips - length(splits)) - 1L)
            * prod(DoubleFactorial(splits + splits - 3L)))
  }
  expect_equal(OldNUS(49:51), NUnrootedSplits(49:51))
  expect_equal(OldNUS(29:32), NUnrootedSplits(29:32))
  expect_equal(OldNUS(1:10), NUnrootedSplits(1:10))
  expect_error(NUnrootedSplits(c(10, 152)))


})

test_that("Log rooted calculated correctly", {
  expect_equal(Log2Rooted(10:30), LnRooted(10:30) / log(2))
  expect_equal(Log2Rooted(10:30), Log2Rooted.int(10:30))
})

test_that("SPR distances calculated correctly", {
  expect_equal(0L, N1Spr(0))
  expect_equal(c(rep(0L, 3), 2L * (3:10 - 3L) * (2L * 3:10 - 7L)),
               vapply(0:10, N1Spr, integer(1)))
  expect_equal(-log2((1L + N1Spr(0:10)) / NUnrooted(0:10)), IC1Spr(0:10))
})
