test_that("PhyToString() works", {
  phy <- StringToPhyDat("012[01]", letters[1:4])
  expect_equal(PhyToString(phy), "012{01}")
  expect_equal(PhyToString(phy, parentheses = "<"), "012<01>")
  expect_equal(PhyToString(phy, parentheses = ">"), "012<01>")
  expect_equal(PhyToString(phy, parentheses = "("), "012(01)")
  expect_equal(PhyToString(phy, parentheses = ")"), "012(01)")
  expect_equal(PhyToString(phy, parentheses = "]"), "012[01]")
  expect_equal(PhyToString(phy, parentheses = "["), "012[01]")
  expect_equal(PhyToString(phy, parentheses = "}"), "012{01}")
  expect_equal(PhyToString(phy, parentheses = "{"), "012{01}")
  expect_equal(PhyToString(phy, parentheses = "!"), "012{01}")
  
  expect_equal(PhyToString(phy, concatenate = FALSE), c(0, 1, 2, "{01}"))
  expect_equal(PhyToString(phy, concatenate = FALSE, ps = ";"),
               paste0(c(0, 1, 2, "{01}"), ";"))
  
  str <- "012{01}0123"
  phy <- StringToPhyDat(str, letters[1:4])
  expect_equal(str, PhyToString(StringToPhyDat(str, letters[1:4])))
  expect_equal(str,
               PhyToString(StringToPhyDat(str, letters[1:4], byTaxon = TRUE),
                           byTaxon = TRUE))
  
  singleChar <- structure(list(a = 1L, b = 1L, c = 2L, d = 2L),
                          weight = 2L, nr = 1L, nc = 2L,
                          index = c(1L, 1L), levels = c("0", "1"),
                          allLevels = c("0", "1", "?"), type = "USER",
                          contrast = structure(c(1, 0, 1, 0, 1, 1), dim = 3:2,
                                               dimnames = list(NULL, c("0", "1"))),
                          class = "phyDat")
  expect_equal(PhyToString(singleChar, concatenate = FALSE, byTaxon = FALSE,
                           useIndex = FALSE), "0011")
})

test_that("PhyToString() supports long levels", {
  ABC <- LETTERS[1:3]
  longLevels <- structure(
    list(
      x = c(8L, 15L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 9L, 10L, 11L, 12L, 13L, 14L),
      y = c(14L, 13L, 12L, 11L, 10L, 9L, 7L, 6L, 5L, 4L, 3L, 2L, 1L, 8L, 15L)),
    weight = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    nr = 15L,
    nc = 14L,
    index = 1:15,
    levels = c(0:6, "-", 7:9, ABC),
    allLevels = c(0:6, "-", 7:9, ABC, "?"),
    type = "USER",
    contrast = structure(c(
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1, 1), dim = 15:14),
    class = "phyDat")
  
  expect_equal(PhyToString(longLevels), "-?0123456789ABCCBA9876543210-?")
  
  expect_equal(PhyToString(longLevels, ps = ";"),
               "-?0123456789ABCCBA9876543210-?;")
  expect_equal(PhyToString(longLevels[, 1], ps = ";", byTaxon = TRUE), "-C;")
  expect_equal(PhyToString(longLevels[1, ], ps = ";", byTaxon = TRUE),
               "-?0123456789ABC;")
  expect_equal(PhyToString(longLevels, ps = ";", useIndex = FALSE,
                           byTaxon = TRUE, concatenate = TRUE),
               "-?0123456789ABCCBA9876543210-?;")
  
  # Two -s → error
  attr(longLevels, "allLevels")[1] <- "-"
  expect_error(PhyToString(longLevels))
})
