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

test_that("UnrootedTreesMatchingSplit works", {
  expect_equal(NRooted(3) * NRooted(5), UnrootedTreesMatchingSplit(c(3, 5)))
  expect_equal(NRooted(30) * NRooted(50), UnrootedTreesMatchingSplit(c(30, 50)))
})

test_that("MultiSplitInformation works", {
  expect_equal(12.8323, MultiSplitInformation(3:5), tolerance = 6)
  expect_equal(CharacterInformation(rep(c('-', '?', 0:2), 1:5)),
               MultiSplitInformation(3:5))
})
