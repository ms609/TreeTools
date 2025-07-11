test_that("PhyToString() supports long levels", {
  skip_if_not_installed("phangorn")
  longLevels <- phangorn::phyDat(rbind(x = c("-", "?", 0:12),
                                       y = c(12:0, "-", "?")),
                                 type = "USER", levels = c(0:6, "-", 7:12))
  expect_equal(PhyToString(longLevels), "-?0123456789ABCCBA9876543210-?")
  
  expect_equal(PhyToString(longLevels, ps = ";"),
               "-?0123456789ABC;CBA9876543210-?;")
  
  # Two -s → error
  attr(longLevels, "allLevels")[1] <- "-"
  expect_error(PhyToString(longLevels))
  
  # 10 → 1
  longLevels <- phangorn::phyDat(rbind(x = c("-", "?", 1:10),
                                       y = c(10:1, "-", "?")),
                                 type = "USER", levels = c(1:6, "-", 7:10))
  expect_equal("-?12345678900987654321-?", PhyToString(longLevels))
})

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
  
  str <- "012{01}0123"
  phy <- StringToPhyDat(str, letters[1:4])
  expect_equal(str, PhyToString(StringToPhyDat(str, letters[1:4])))
  expect_equal(str,
               PhyToString(StringToPhyDat(str, letters[1:4], byTaxon = TRUE),
                           byTaxon = TRUE))
})
