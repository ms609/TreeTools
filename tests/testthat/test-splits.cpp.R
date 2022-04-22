test_that("Bad data results in error", {
  skip_if(Sys.getenv("USING_ASAN") != "")
  expect_error(cpp_edge_to_splits(matrix(1, 3, 3), 1:3, 3),
               "must contain two col")
  expect_error(cpp_edge_to_splits(matrix(1, 10, 2), 1:10, 0),
               "must contain tips")
  expect_error(cpp_edge_to_splits(matrix(1, 10, 2), 1:10, -10),
               "must contain tips")
  expect_error(cpp_edge_to_splits(matrix(1, 10, 2), 1:9, 5),
               "ength of `order` must equal number of edges")
  expect_error(cpp_edge_to_splits(matrix(1, 2, 2), 1:2, 6),
               "Not enough edges") # Impossible anyway...
  
  s9 <- as.Splits(BalancedTree(9))
  expect_error(duplicated_splits(structure(s9, nTip = NULL), TRUE),
               "`nTip` attribute")
  expect_error(duplicated_splits(structure(s9, nTip = 999L), TRUE),
               "tip number")
  
  expect_error(s9 & s9[[1]], "same number of splits")
  expect_error(s9 | s9[[1]], "same number of splits")
  expect_error(xor(s9, s9[[1]]), "same number of splits")
  
  noNTip <- structure(s9, nTip = NULL)
  expect_error(noNTip & s9, "`x` lacks nTip attrib")
  expect_error(noNTip | s9, "`x` lacks nTip attrib")
  expect_error(xor(noNTip, s9), "`x` lacks nTip attrib")
  
  expect_error(s9 & noNTip, "`y` lacks nTip attrib")
  expect_error(s9 | noNTip, "`y` lacks nTip attrib")
  expect_error(xor(s9, noNTip), "`y` lacks nTip attrib")

  wrongTip <- structure(s9, nTip = 15L)
  expect_error(xor(s9, wrongTip), "`y` differ in `nTip`")
  expect_error(s9 & wrongTip, "`y` differ in `nTip`")
  expect_error(s9 | wrongTip, "`y` differ in `nTip`")
})
