nasty <- structure(list(edge = structure(
  c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
    5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
  .Dim = c(12, 2)),
  Nnode = 5L,
  tip.label = letters[1:8]),
  class = "phylo") # Danger: Do not plot!

test_that("keep_tip() fails with small edge", {
  # False positives with ASaN
  skip_if(Sys.getenv("USING_ASAN") != "")
  expect_error(keep_tip(BalancedTree(8)$edge[, 1, drop = FALSE],
                        !tabulate(6:4, 8)))
  expect_error(keep_tip(BalancedTree(8)$edge[, 1], !tabulate(6:4, 8)))
})

test_that("keep_tip() works", {
  expect_error(keep_tip(BalancedTree(8)$edge[, c(1, 2, 1)], !tabulate(6:4, 8)))
  
  expect_equal(keep_tip(BalancedTree(9)$edge, !tabulate(5:8, 9)),
               matrix(c(6, 7, 8, 9, 9, 8, 7, 6,
                        7, 8, 9, 1, 2, 3, 4, 5), 8, 2))
  
  expect_equal(keep_tip(BalancedTree(8)$edge, !tabulate(6:4, 8)),
               matrix(c(6, 7, 8, 8, 7, 6, 9, 9,
                        7, 8, 1, 2, 3, 9, 4, 5), 8, 2))
  expect_equal(keep_tip(BalancedTree(8)$edge, !tabulate(1:4, 8)),
               BalancedTree(4)$edge)
  
  expect_equal(keep_tip(BalancedTree(8)$edge, !tabulate(3:8, 8)),
               BalancedTree(2)$edge)
  
  expect_equal(keep_tip(UnrootTree(BalancedTree(4))$edge, !tabulate(1, 4)),
               matrix(c(4, 4, 4, 1, 2, 3), 3, 2))
  
  expect_equal(Preorder(keep_tip(as.phylo(3, 7)$edge, !tabulate(1:4, 7))),
               PectinateTree(3)$edge)
  
  # Rooted tree may lose root when reaching a polytomy:
  expect_equal(keep_tip(root(StarTree(6), 4, resolve.root = TRUE)$edge,
                        !tabulate(4, 6)),
               StarTree(5)$edge)
  
  # But need not:
  expect_equal(keep_tip(RootTree(BalancedTree(5), 3)$edge, !tabulate(3, 5)),
               BalancedTree(4)$edge)
  
  unrooted <- ape::read.tree(text = "(a, b, (c, d, ((e1, e2), (f, g))));")
  expect_equal(keep_tip(unrooted$edge, !tabulate(1:4, 8)),
               UnrootTree(BalancedTree(4))$edge)
})

test_that("DropTip() works", {
  bal8 <- BalancedTree(8)
  expect_equal(NTip(DropTip(bal8, 1:8, check = FALSE)), 0)
  expect_warning(expect_true(all.equal(bal8, DropTip(bal8, -1))),
                 "`tip` must be > 0")
  expect_warning(expect_true(all.equal(bal8, DropTip(bal8, 99))),
                 "Tree only has 15 nodes")
  expect_warning(expect_true(all.equal(bal8, DropTip(bal8, "MissingTip"))),
                 "not present in tree")
  expect_warning(expect_identical(
    DropTip(bal8, c("t8", "NotThere"), check = TRUE),
    DropTip(bal8, c("t8", "NotThere"), check = FALSE)), "not present in tree")
  expect_error(DropTip(bal8, TRUE),
               "`tip` must list `TRUE` or `FALSE` for each leaf")
  expect_error(DropTip(bal8, list("Invalid format")), "`tip` must be of type")
  
  expect_equal(DropTip(bal8, 7:8), DropTip(bal8, 15L))
  expect_true(all.equal(ape::drop.tip(bal8, 6:8), DropTip(bal8, 6:8)))
  expect_true(all.equal(ape::drop.tip(bal8, c(3, 5, 7)), DropTip(bal8, c(3, 5, 7))))
  
  expect_equal(Preorder(DropTip(Preorder(nasty), c(1, 3))),
               Preorder(DropTip(nasty, c(1, 3))))
})

test_that("DropTip() root relocation", {
  nTip <- 12
  nKept <- nTip / 2L
  for (i in c(8, 30, 324)) { # A selection of interesting failure cases
    set.seed(i)
    
    bigTree <- RandomTree(nTip)
    bigDrop <- sample.int(nTip, nKept)
    bigKeep <- setdiff(seq_len(nTip), bigDrop)
    
    reduced <- DropTip(bigTree, bigDrop)
    expect_true(all.equal(reduced, unroot(drop.tip(bigTree, bigDrop))))
    
    map <- which(KeptVerts(bigTree, !tabulate(bigDrop, nTip)))
    expect_equal(map[seq_len(nKept)], sort(bigKeep))
    newDesc <- .ListDescendants(reduced)
    bigDesc <- .ListDescendants(bigTree)
    
    for (newNode in (nKept + 1):(nKept * 2 - 2)) {
      expect_equal(sort(map[newDesc[[newNode]]]),
                   sort(intersect(bigDesc[[map[newNode]]], bigKeep)))
    }
  }
  #microbenchmark::microbenchmark(ape::drop.tip(bigTree, bigTip), DropTip(bigTree, bigTip))
  
  nTip <- 1284
  set.seed(43)
  bigTree <- RandomTree(nTip)
  bigDrop <- sample.int(nTip, nKept)
  bigKeep <- setdiff(seq_len(nTip), bigDrop)
  
  reduced <- DropTip(bigTree, bigDrop)
  expect_true(all.equal(reduced, unroot(drop.tip(bigTree, bigDrop))))
})

test_that("DropTip.multiPhylo() with attributes", {
  multi <- c(bal8 = BalancedTree(8), pec8 = PectinateTree(8))
  attr(multi, "TipLabel") <- paste0("t", 1:8)
  
  expect_equal(DropTip(unclass(multi), "t6"), DropTip(multi, "t6"))
  
  expect_equal(attr(DropTip(multi, "t8"), "TipLabel"),
               paste0("t", 1:7))
  expect_equal(names(DropTip(multi, "t8")), names(multi))
  expect_equal(DropTip(multi[1], "t1")[[1]], DropTip(multi[[1]], "t1"))
})

test_that("DropTip.Splits()", {
  bal9 <- BalancedTree(9)
  s9 <- as.Splits(bal9)
  s19 <- as.Splits(PectinateTree(19))
  
  expect_error(DropTip(s9, c(T, F)), "each leaf\\b")
  expect_warning(
    expect_warning(
      expect_equal(DropTip(s9, c(0, 9:10)), DropTip(s9, !tabulate(1:8, 9))),
      "only has 9 leaves"),
    "`tip` must be > 0"
  )
  expect_error(DropTip(s9, raw(1)), "`tip` must be of type")
  
  expect_equal(s9, DropTip(s9, NULL))
  expect_warning(
    expect_equal(DropTip(s9, c("t7", "missing")),
                 DropTip(s9, "t7", check = FALSE)),
    "not present in tree")
  
  expect_equal(unname(DropTip(s9, 4:5)), unname(as.Splits(DropTip(bal9, 4:5))))
  expect_equal(unname(KeepTip(s9, c(1:4, 7:9))),
               unname(as.Splits(DropTip(bal9, 6:5))))
  
  expect_equal(DropTip(s9[[1:5]], 8:9),
               DropTip(s9, 8:9))
  
  expect_equal(thin_splits(s9, !logical(9)),
                           structure(raw(0), .Dim = c(0L, 0L)))
  expect_equal(thin_splits(s19, 1:19 %in% 2:19),
               structure(raw(0), .Dim = c(0L, 1L)))
  expect_equal(thin_splits(s9, logical(9)), s9, ignore_attr = TRUE)
})

test_that("KeepTip() works", {
  expect_warning(expect_true(all.equal(
    BalancedTree(paste0("t", 5:8)),
    KeepTip(BalancedTree(8), paste0("t", 5:9))
  )))
  
  expect_warning(expect_true(all.equal(
    BalancedTree(paste0("t", 5:8)),
    KeepTip(BalancedTree(8), 5:9)
  )))
  
  expect_true(all.equal(
    KeepTip(Postorder(UnrootTree(BalancedTree(letters[1:12]))), 6:12),
    read.tree(text = "(f, (i, (g, h)), (l, (j, k)));")
  ))
  
  expect_identical(
    KeepTip(BalancedTree(12), !tabulate(1:5, 12)),
    DropTip(BalancedTree(12), !tabulate(6:12, 12))
    )
})

test_that("KeepTip() retains edge lengths", {
  bal9 <- BalancedTree(9)
  bal9$edge.length <- 10 * 1:16
  expect_equal(KeepTip(bal9, c(1, 3, 5, 7, 9))[["edge.length"]],
               c(10, 20, 30 + 40,
                 60,
                 70 + 90,
                 100, 110 + 130,
                 140 + 160))
})

test_that("KeepTipPreorder()/Postorder()", {
  pre <- Preorder(BalancedTree(1:9))
  post <- Postorder(pre)
  keep <- !tabulate(7, 9)
  expect_true(all.equal(KeepTipPreorder(pre, keep),
                        KeepTip(pre, keep)))
  expect_true(all.equal(KeepTipPostorder(post, keep),
                        KeepTip(post, keep)))
  
  emptyTree <- structure(list(edge = matrix(0, 0, 2),
                              tip.label = character(0),
                              Nnode = 0), class = "phylo")
  expect_equal(KeepTipPreorder(pre, logical(9)),
               emptyTree)
  expect_equal(KeepTipPostorder(post, logical(9)),
               emptyTree)

  expect_equal(KeepTipPreorder(pre, !logical(9)),
               pre)
  expect_equal(KeepTipPostorder(post, !logical(9)),
               post)
})
