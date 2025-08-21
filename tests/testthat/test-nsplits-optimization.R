test_that("NSplits optimization produces same results", {
  # Test trees with different structures
  
  # Helper function: Original NSplits implementation using collapse.singles concept
  NSplits_original <- function(x) {
    if (length(x[["tip.label"]]) < 4L) {
      0L
    } else {
      # TreeIsRooted implementation
      tree_is_rooted <- function(tree) {
        edge <- tree[["edge"]]
        parent <- edge[, 1]
        sum(parent == min(parent)) < 3L
      }
      
      # For binary trees without singles, Nnode should equal the number of internal nodes
      # This simulates what collapse.singles would return for a tree without single nodes
      x[["Nnode"]] - 1L - tree_is_rooted(x)
    }
  }
  
  # Test case 1: Simple 4-tip binary tree ((t1,t2),(t3,t4))
  tree4 <- structure(list(
    edge = matrix(c(5, 1, 5, 2, 6, 3, 6, 4, 7, 5, 7, 6), ncol = 2, byrow = TRUE),
    tip.label = c("t1", "t2", "t3", "t4"),
    Nnode = 3L
  ), class = "phylo")
  
  # Test case 2: 6-tip binary tree (((t1,t2),(t3,t4)),(t5,t6))
  tree6 <- structure(list(
    edge = matrix(c(7, 1, 7, 2, 8, 3, 8, 4, 9, 5, 9, 6, 10, 7, 10, 8, 11, 9, 11, 10), 
                  ncol = 2, byrow = TRUE),
    tip.label = c("t1", "t2", "t3", "t4", "t5", "t6"),
    Nnode = 5L
  ), class = "phylo")
  
  # Test case 3: Small trees (< 4 tips) should return 0
  tree2 <- structure(list(
    edge = matrix(c(3, 1, 3, 2), ncol = 2, byrow = TRUE),
    tip.label = c("t1", "t2"),
    Nnode = 1L
  ), class = "phylo")
  
  tree3 <- structure(list(
    edge = matrix(c(4, 1, 4, 2, 4, 3), ncol = 2, byrow = TRUE),
    tip.label = c("t1", "t2", "t3"),
    Nnode = 1L
  ), class = "phylo")
  
  # Test case 4: Larger tree (8 tips) - pectinate tree
  tree8 <- structure(list(
    edge = matrix(c(9, 1, 10, 2, 11, 3, 12, 4, 13, 5, 14, 6, 15, 7, 15, 8,
                    14, 15, 13, 14, 12, 13, 11, 12, 10, 11, 9, 10), 
                  ncol = 2, byrow = TRUE),
    tip.label = paste0("t", 1:8),
    Nnode = 7L
  ), class = "phylo")
  
  # Test our optimized version against the original logic
  expect_equal(NSplits(tree4), NSplits_original(tree4))
  expect_equal(NSplits(tree6), NSplits_original(tree6))
  expect_equal(NSplits(tree8), NSplits_original(tree8))
  
  # Small trees should return 0
  expect_equal(NSplits(tree2), 0L)
  expect_equal(NSplits(tree3), 0L)
  
  # Specific expected values based on the formula: Nnode - 1 - TreeIsRooted
  # tree4: 3 - 1 - 1 = 1 (rooted)
  # tree6: 5 - 1 - 1 = 3 (rooted) 
  # tree8: 7 - 1 - 1 = 5 (rooted)
  expect_equal(NSplits(tree4), 1L)
  expect_equal(NSplits(tree6), 3L)
  expect_equal(NSplits(tree8), 5L)
})

test_that("cpp_count_splits handles edge cases correctly", {
  # Test the C++ function directly
  
  # Empty edge matrix
  empty_edge <- matrix(integer(0), ncol = 2)
  expect_equal(cpp_count_splits(empty_edge, 0L), 0L)
  
  # 2 tips tree
  edge_2tips <- matrix(c(3, 1, 3, 2), ncol = 2, byrow = TRUE)
  expect_equal(cpp_count_splits(edge_2tips, 2L), 0L)
  
  # 3 tips tree
  edge_3tips <- matrix(c(4, 1, 4, 2, 4, 3), ncol = 2, byrow = TRUE)
  expect_equal(cpp_count_splits(edge_3tips, 3L), 0L)
  
  # Invalid edge matrix (wrong number of columns)
  invalid_edge <- matrix(c(1, 2, 3), ncol = 3)
  expect_error(cpp_count_splits(invalid_edge, 2L), "Edge matrix must contain two columns")
})