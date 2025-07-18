test_that("Reweight()", {
  mat <- rbind(
    a = c(c1 = 0, c2 = 0, c3 = 2, c4 = 0),
    b = c(0, 0, 2, 0),
    c = c(1, 0, 3, 0),
    d = c(1, 0, 3, 0))
  dat <- MatrixToPhyDat(mat)
  
  expect_error(Reweight(mat, letters[1:4]),
               "must be a numeric vector")
  expect_error(Reweight(mat, 1:2),
               "Length of `weights` \\(2\\) must match .* characters \\(4\\)")
  expect_error(Reweight(mat, 1:8),
               "Length of `weights` \\(8\\) must match .* characters \\(4\\)")
  expect_error(Reweight(mat, c("-1" = 2, "999" = 911)),
               "Indices in `names\\(weights\\)` .*ween 1 and 4. Found -1, 999")
  
  newWeights <- c(0, 1, 3, 1)
  expect_equal(
    Reweight(mat, newWeights),
    `colnames<-`(mat[, c(2, 3, 3, 3, 4)],
                 c("c2_1", "c3_1", "c3_2", "c3_3", "c4_1"))
  )
  
  rwDat <- Reweight(dat, c("3" = 3, "2" = 1, "1" = 0))
  expect_equal(PhyDatToMatrix(rwDat), PhyDatToMatrix(dat))
  rwAtt <- attributes(rwDat)
  sameAtt <- c("names", "nc", "levels", "allLevels", "contrast", "type",
               "nr", "index", "class")
  expect_equal(rwAtt[sameAtt], attributes(dat)[sameAtt])
  expect_equal(sum(rwAtt[["weight"]]), sum(newWeights))
  expect_equal(rwAtt[["weight"]], c(0, 2, 3))
})

test_that("Reweight() okay with tree search", {
  skip_if_not_installed("TreeSearch")
  library("TreeSearch")
  tree <- NJTree(Lobo.phy)
  expect_equal(TreeLength(tree, Reweight(Lobo.phy,
                                       rep_len(2, length(Lobo.data[[1]])))),
               TreeLength(tree, Lobo.phy) * 2)
  
  expect_equal(TreeLength(tree, Reweight(Lobo.phy[, 9:12], c(0, 100, 1, 0))),
               TreeLength(tree, Lobo.phy[, 10]) * 100 +
                 TreeLength(tree, Lobo.phy[, 11])) # = 105
})
