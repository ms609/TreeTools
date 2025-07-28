test_that("Decompose()", {
  data("Lobo")
  
  expect_error(Decompose("Lobo.phy"), "phyDat")

  expect_error(Decompose(Lobo.phy, logical(11)),
               "length.indices.* number of characters in .dataset.")
  
  expect_warning(decomposed <- Decompose(Lobo.phy, c(1234:1235, 11)), 
                 "115 characters; .*1234, 1235 not found")

  NumberOfChars <- function(x) sum(attr(x, "weight"))
  expect_equal(NumberOfChars(decomposed),
               NumberOfChars(Lobo.phy) + 1)
  
  expect_equal(attr(decomposed, "originalIndex"),
               c(1:11, 11:115))
               
  expect_equal(attr(Decompose(Lobo.phy, 10:12), "originalIndex"),
               c(1:11, 11:115))
  
  decompMat <- as.matrix(decomposed)
  taxa <- c(8, 11, 12, 37, 36)
  if (interactive()) {
    dput(as.character(as.matrix(Lobo.phy)[taxa, 11]))
    # = c("?", "{01}", "0", "1", "2")
  }
  expect_equal(decompMat[taxa, 11:12],
               matrix(c("?", "{01}", "0", "1", "1",
                        "?", "0", "0", "0", "1"), ncol = 2,
                      dimnames = list(names(Lobo.phy)[taxa])))
  
  unchanged <- Decompose(Lobo.phy, logical(115))
  
  expect_equal(gsub("\\{\\-(\\d)\\}", "{\\1-}", perl = TRUE,
                    as.matrix(unchanged)), # Order unimportant
               as.matrix(Lobo.phy))
  
  expect_equal(Decompose(Lobo.phy[, 1], TRUE),
               structure(Lobo.phy[, 1], originalIndex = 1))
  
})

test_that("Decompose() with ambiguities", {
  gme <- structure(list(
    tax1 = c(1L, 4L, 7L), tax2 = c(2L, 5L, 4L),
    tax3 = c(3L, 6L, 8L), tax4 = c(2L, 6L, 9L)),
    weight = c(1L, 1L, 1L), nr = 3L, nc = 10L, index = 1:3,
    levels = c("-", "0", "1", "2", "3", "4", "5", "6", "7", "8"),
    allLevels = c("{0123}", "{24}", "2", "0", "{5678}", "5", "-", "{01}", "7"),
    type = "USER",
    contrast = structure(
      c(0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0,
        0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0), dim = 9:10, 
      dimnames = list(NULL, c("-", 0:8))),
    class = "phyDat")
  
  # "{0123}{24}2{24}"
  exp1 <- c("???0", "11??", "1100", "11??")
  expect_equal(PhyToString(Decompose(gme[, 1], 1), concatenate = FALSE), exp1)
  
  # "0{5678}55"
  exp2 <- c("00000000", "11111{01}{01}{01}", "11111000", "11111000")
  # Need to include char 1 so levels aren't lost from contrast when subsetting
  expect_equal(PhyToString(Decompose(gme[, 1:2], 2), concatenate = FALSE),
               paste0(PhyToString(gme[, 1], concat = FALSE), exp2))
  
  # "-0{01}7"
  exp3 <- c("-------", "0000000", "{01}000000", "1111111")
  expect_equal(PhyToString(Decompose(gme, 3), concatenate = FALSE),
               paste0(PhyToString(gme[, 1:2], concatenate = FALSE), exp3))
  expect_equal(PhyToString(Decompose(gme, 1:3), concatenate = FALSE),
               gsub("?", "{01}", fixed = TRUE, paste0(exp1, exp2, exp3)))
})
