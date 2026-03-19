test_that("cpp_edge_to_splits_batch() matches per-tree results", {
  trees <- as.phylo(c(30899669, 9149275, 12823175, 19740197), 11L,
                    seq_len(11L))
  nTip <- 11L
  edges <- lapply(trees, `[[`, "edge")
  orders <- lapply(trees, function(tr) {
    postorder_order(tr[["edge"]]) - 1L
  })
  
  batchResult <- cpp_edge_to_splits_batch(edges, orders, nTip)
  perTree <- lapply(seq_along(trees), function(i) {
    cpp_edge_to_splits(edges[[i]], orders[[i]], nTip)
  })
  
  expect_length(batchResult, length(trees))
  for (i in seq_along(trees)) {
    expect_equal(batchResult[[i]], perTree[[i]])
  }
})

test_that("cpp_edge_to_splits_batch() handles edge cases", {
  # Empty list
  expect_length(cpp_edge_to_splits_batch(list(), list(), 5L), 0)

  # n_tip == 0
  result <- cpp_edge_to_splits_batch(
    list(matrix(1L, 10, 2)), list(0:9), 0L
  )
  expect_equal(result[[1]], matrix(raw(0), 0, 0))
  
  # Negative n_tip
  expect_error(cpp_edge_to_splits_batch(list(), list(), -1L),
               "non-negative number of tips")
  
  # Mismatched list lengths
  expect_error(
    cpp_edge_to_splits_batch(list(matrix(1L, 3, 2)), list(), 5L),
    "same length"
  )
  
  # Bad edge matrix (3 columns)
  expect_error(
    cpp_edge_to_splits_batch(list(matrix(1L, 3, 3)), list(0:2), 3L),
    "two columns"
  )
  
  # Mismatched order length
  expect_error(
    cpp_edge_to_splits_batch(list(matrix(1L, 10, 2)), list(0:8), 5L),
    "number of edges"
  )
})

test_that("Bad data results in error", {
  skip_if(Sys.getenv("USING_ASAN") != "")
  expect_error(cpp_edge_to_splits(matrix(1, 3, 3), 1:3, 3),
               "must contain two col")
  expect_equal(cpp_edge_to_splits(matrix(1, 10, 2), 1:10, 0),
               matrix(raw(0), 0, 0))
  expect_error(cpp_edge_to_splits(matrix(1, 10, 2), 1:10, -10),
               "Tree must contain non-negative number of tips")
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
