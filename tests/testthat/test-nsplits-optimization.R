test_that("NSplits optimization produces same results", {
  # Test trees with different structures
  
  # Original NSplits implementation
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
      
      ape::collapse.singles(x)[["Nnode"]] - 1L - TreeIsRooted(x)
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
  
  # Test case 5: Tree with singleton nodes (degree-2 nodes)
  # Represents tree equivalent to (a, (b, (c), (d, e))) after parsing
  # but with an extra internal node (singleton) connecting to c
  tree_with_singles <- structure(list(
    edge = matrix(c(
      6, 1,   # root -> a
      6, 7,   # root -> internal
      7, 2,   # internal -> b  
      7, 8,   # internal -> singleton node
      8, 3,   # singleton -> c (this makes node 8 a singleton)
      7, 9,   # internal -> another internal
      9, 4,   # internal -> d
      9, 5    # internal -> e
    ), ncol = 2, byrow = TRUE),
    tip.label = c("a", "b", "c", "d", "e"),
    Nnode = 4L  # 4 internal nodes total: 6 (root), 7, 8 (singleton), 9
  ), class = "phylo")
  
  # Test our optimized version against the original logic for binary trees
  expect_equal(NSplits(tree4), NSplits_original(tree4))
  expect_equal(NSplits(tree6), NSplits_original(tree6))
  expect_equal(NSplits(tree8), NSplits_original(tree8))
  
  # Small trees should return 0
  expect_equal(NSplits(tree2), 0L)
  expect_equal(NSplits(tree3), 0L)
  
  # Tree with singleton nodes should have fewer splits after collapsing singletons
  # Expected: node 8 is singleton (appears once as parent [8->3], once as child [7->8])
  # After collapsing: we have 3 remaining internal nodes (6, 7, 9)
  # Formula: (3 internal - 0 remaining singleton) - 1 - 1(rooted) = 1 split
  expect_equal(NSplits(tree_with_singles), 1L)
  
  # Specific expected values based on the formula: (Nnode - n_singles) - 1 - TreeIsRooted
  # tree4: (3 - 0) - 1 - 1 = 1 (rooted, no singletons)
  # tree6: (5 - 0) - 1 - 1 = 3 (rooted, no singletons) 
  # tree8: (7 - 0) - 1 - 1 = 5 (rooted, no singletons)
  expect_equal(NSplits(tree4), 1L)
  expect_equal(NSplits(tree6), 3L)
  expect_equal(NSplits(tree8), 5L)
})

test_that("NSplits handles trees with non-preorder edge numbering", {
  # Test tree with "nasty" edge order (not preorder)
  nasty_edge <- structure(c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
                           5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
                         .Dim = c(12, 2))
  nasty_tree <- structure(list(edge = nasty_edge, Nnode = 5L, tip.label = letters[1:8]),
                         class = "phylo")
  
  expect_equal(NSplits(nasty_tree), 3L)
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
  
  # Test tree with singleton: 4 tips but with singleton node
  # Tree: (a, (b, (c), d)) - node connecting just c is singleton
  edge_singleton <- matrix(c(
    5, 1,   # root -> a
    5, 6,   # root -> internal 
    6, 2,   # internal -> b
    6, 7,   # internal -> singleton
    7, 3,   # singleton -> c (node 7 is singleton)
    6, 4    # internal -> d
  ), ncol = 2, byrow = TRUE)
  
  # After collapsing singleton node 7: we get (a, (b, c, d))
  # This should have 2 internal nodes (5=root, 6=internal)
  # Result: (2 - 0) - 1 - 1 = 0 splits
  expect_equal(cpp_count_splits(edge_singleton, 4L), 0L)
  
  # Invalid edge matrix (wrong number of columns)
  invalid_edge <- matrix(c(1, 2, 3), ncol = 3)
  expect_error(cpp_count_splits(invalid_edge, 2L), "Edge matrix must contain two columns")
})
