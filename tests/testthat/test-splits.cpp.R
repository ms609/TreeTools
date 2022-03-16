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
  
  expect_error(xor(splits, splits[[1]]),
               "same number of splits")
  expect_error(xor(structure(splits, nTip = NULL), splits),
               "`x` lacks nTip attrib")
  expect_error(xor(splits, structure(splits, nTip = 15)),
               "`y` differ in `nTip`")
})
