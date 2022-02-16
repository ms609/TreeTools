test_that("KeptVerts() works", {
  Keepers <- function (n, nTip) as.logical(tabulate(n, nTip))
  
  bal12 <- BalancedTree(12)
  Y <- TRUE
  N <- FALSE
  expect_equal(KeptVerts(bal12, 1:12 %in% 1:9),
               c(1:12 %in% 1:9, rep(Y, 6), N, Y, Y, N, N))
  expect_equal(KeptVerts(bal12, 1:12 %in% 7:12),
               c(1:12 %in% 7:12, rep(N, 6), rep(Y, 5)))
  expect_equal(KeptVerts(bal12, 1:12 %in% 4:6),
               c(1:12 %in% 4:6, rep(N, 4), rep(Y, 2), rep(N, 5)))
  expect_equal(KeptVerts(bal12, 1:12 %in% c(1, 2, 6, 12)),
               c(1:12 %in% c(1, 2, 6, 12), Y, Y, N, Y, N, N, N, N, N, N, N))
  
  pec9 <- UnrootTree(PectinateTree(9))
  kept <- c(1, 3, 8, 9)
  expect_equal(which(KeptVerts(pec9, 1:9 %in% kept)),
               c(kept, 10, 16))
  
  real <- read.tree(text = 
  "(t1:3,t2:1,((t3:1,t4:1):3,(t5:2,(t6:1,(t7:6,((t8:8,(t9:2,t10:4):3):2,(t11:1,t12:2):5):1):2):1):1):2);")
  kept <- c(1L, 3L, 4L, 8L, 9L)
  # Duplicate root is not retained.
  expect_equal(
    which(KeptVerts(real, Keepers(kept, 12))),
    c(kept, 13, 15, 20))
  
  # Duplicate keep_tip tests
  expect_equal(which(KeptVerts(BalancedTree(9)$edge, !tabulate(5:8, 9))),
               c(1:4, 9, 10:13))
  
  expect_equal(which(KeptVerts(BalancedTree(8)$edge, !tabulate(6:4, 8))),
               c(1:3, 7:8, 9:11, 15))
  expect_equal(which(KeptVerts(BalancedTree(8)$edge, !tabulate(1:4, 8))),
               c(5:8, 13:15))
  
  expect_equal(which(KeptVerts(BalancedTree(8)$edge, !tabulate(3:8, 8))),
               c(1:2, 11))
  
  expect_equal(which(KeptVerts(ape::unroot(BalancedTree(4)), !tabulate(1, 4))),
               c(2:4, 5))
  
  expect_equal(which(KeptVerts(as.phylo(3, 7)$edge, !tabulate(1:4, 7))),
               c(5:7, 9, 11))
  
  # Rooted tree may lose root when reaching a polytomy:
  expect_equal(which(KeptVerts(root(StarTree(6), 4, resolve.root = TRUE),
                               !tabulate(4, 6))),
               c(1:3, 5:6, 8))
  
  # But need not:
  expect_equal(which(KeptVerts(RootTree(BalancedTree(5), 3)$edge,
                               !tabulate(3, 5))),
               c(1:2, 4:5, 7:9))
  
})
