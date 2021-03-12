context("Information.R")

test_that("Trees matching splits calculated correctly", {
  expect_equal(NUnrooted(9), TreesMatchingSplit(0, 9))
  expect_equal(NUnrooted(9), TreesMatchingSplit(9, 0))
  expect_equal(NUnrooted(9), TreesMatchingSplit(8, 1))
  expect_equal(LnUnrooted(9), LnTreesMatchingSplit(0, 9))
  expect_equal(LnUnrooted.int(9), LnTreesMatchingSplit(9, 0))
  expect_equal(LnUnrooted.int(9), LnTreesMatchingSplit(1, 8))
  expect_equal(NUnrooted(9), TreesMatchingSplit(8, 1))
  expect_equal(Log2Unrooted(9), Log2TreesMatchingSplit(0, 9))
  expect_equal(Log2Unrooted.int(9), Log2TreesMatchingSplit(9, 0))
  expect_equal(Log2Unrooted.int(9), Log2TreesMatchingSplit(1, 8))
  expect_equal(log(315/10395)/-log(2), SplitInformation(3, 5))
})

test_that("UnrootedTreesMatchingSplit() correct", {
  expect_equal(NRooted(3) * NRooted(5), UnrootedTreesMatchingSplit(c(3, 5)))
  expect_equal(LnRooted(30) + LnRooted(50), LnUnrootedTreesMatchingSplit(30, 50))
  expect_equal(sum(log2DoubleFactorials[15],
                   Log2DoubleFactorial(1:4 + 1:4 - 3L))
               - log2DoubleFactorials[11],
               Log2UnrootedTreesMatchingSplit(1:4))
})

test_that("MultiSplitInformation() works", {
  expect_equal(12.8323, MultiSplitInformation(3:5), tolerance = 6)
  expect_equal(CharacterInformation(rep(c('-', '?', 0:2), 1:5)),
               MultiSplitInformation(3:5))
})

test_that("TreesMatchingSplit() accepts different formats", {
  expect_equal(TreesMatchingSplit(4, 5), TreesMatchingSplit(4:5))
  expect_equal(LnTreesMatchingSplit(4, 5), LnTreesMatchingSplit(4:5))
  expect_equal(Log2TreesMatchingSplit(4, 5), Log2TreesMatchingSplit(4:5))
})
