test_that("Decompose()", {
  data("Lobo")

  expect_error(Decompose(Lobo.phy, logical(11)),
               "length.indices.* number of characters in .dataset.")
  
  expect_warning(decomposed <- Decompose(Lobo.phy, c(1234:1235, 11)), 
                 "115 characters; .*1234, 1235 not found")

  NumberOfChars <- function(x) sum(attr(x, "weight"))
  expect_equal(NumberOfChars(decomposed),
               NumberOfChars(Lobo.phy) + 1)
  
  decompMat <- as.matrix(decomposed)
  taxa <- c(8, 11, 12, 37, 36)
  if (interactive()) {
    dput(as.character(as.matrix(Lobo.phy)[taxa, 11]))
    # = c("?", "{01}", "0", "1", "2")
  }
  expect_equal(decompMat[taxa, 11:12],
               matrix(c("?", "{01}", "0", "1", "1",
                        "?", "{0}", "0", "0", "1"), ncol = 2,
                      dimnames = list(names(Lobo.phy)[taxa])))
  
  expect_equal(Decompose(Lobo.phy, logical(115)), Lobo.phy)
  
  expect_equal(Decompose(Lobo.phy[, 1], TRUE), Lobo.phy[, 1])
  
})
