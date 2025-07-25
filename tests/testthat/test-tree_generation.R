test_that("Pectinate trees are generated", {
  expect_equal(ape::read.tree(text = "(t1, (t2, (t3, t4)));"),
               PectinateTree(4L))
  expect_equal(ape::read.tree(text = "(a, (b, (c, (d, e))));"),
               PectinateTree(letters[1:5]))
  expect_equal(ape::read.tree(text = "(a, (b, (c, (d, e))));"),
               PectinateTree(ape::read.tree(text = "(a, ((b, c), (d, e)));")))
  data("Lobo", package = "TreeTools")
  expect_equal(ape::read.tree(text = "(Cricocosmia, (Aysheaia, Siberion));"),
               PectinateTree(.SubsetPhyDat(Lobo.phy, 2:4)))
  expect_true(is.integer(PectinateTree(8)$edge))
})

test_that("Balanced trees are generated correctly", {
  expect_equal(BalancedTree(0), ZeroTaxonTree())
  # nTip even
  expect_equal(BalancedTree(8L),
               Preorder(ape::read.tree(
                 text = "(((t1, t2), (t3, t4)), ((t5, t6), (t7, t8)));")
               ))
  # nTip odd
  expect_equal(BalancedTree(9L),
               Preorder(ape::read.tree(
                 text = "((((t1, t2), t3), (t4, t5)), ((t6, t7), (t8, t9)));")
               ))
  expect_equal(BalancedTree(9L)$edge, Preorder(BalancedTree(9L)$edge))
  expect_equal(BalancedTree(as.character(1:9)), BalancedTree(1:9))
  escapees <- c("Apostrophe's", 'and quote"s')
  expect_equal(ignore_attr = TRUE,
               PectinateTree(escapees), BalancedTree(escapees))
  expect_equal(integer(0), .BalancedBit(seq_len(0)))
  expect_equal("Test", .BalancedBit("Test"))
  expect_true(is.integer(BalancedTree(8)[["edge"]]))
})

test_that("BalancedTree(lengths)", {
  expect_equal(BalancedTree(0, lengths = 1), ZeroTaxonTree())
  expect_equal(BalancedTree(1, lengths = 2), SingleTaxonTree("t1", lengths = 2))
  expect_equal(BalancedTree(2, lengths = 2)$edge.length, c(2, 2))
  expect_equal(BalancedTree(3, lengths = 1:3)$edge.length, c(1:3, 1))
})
test_that("PectinateTree(lengths)", {
  expect_equal(PectinateTree(0, lengths = 1), ZeroTaxonTree())
  expect_equal(PectinateTree(1, lengths = 2), SingleTaxonTree("t1", lengths = 2))
  expect_equal(PectinateTree(2, lengths = 2)$edge.length, c(2, 2))
  expect_equal(PectinateTree(3, lengths = 1:3)$edge.length, c(1:3, 1))
})

test_that("StarTree() works", {
  expect_equal(ape::read.tree(text = "(t1, t2, t3, t4, t5, t6, t7, t8);"),
               StarTree(8L))
  expect_true(is.integer(StarTree(8)$edge))
  expect_null(StarTree(8L)[["edge.length"]])
  expect_equal(StarTree(8L, 8:1)[["edge.length"]], 8:1)
})

test_that("Random trees are generated correctly", {
  # These two seeds give two different node labellings
  random3 <- c(4, 4, 5, 5, 1, 5, 2, 3)
  set.seed(1)
  expect_equal(random3, RandomTree(3, root = 1)$edge[1:8])
  set.seed(4)
  expect_equal(random3, RandomTree(3, root = "t1")$edge[1:8])
  
  expect_true(all.equal(RandomTree(3, root = "t2"),
    PectinateTree(c("t2", "t3", "t1"))))
  expect_equal(c(4, 4, 4), RandomTree(3, root = FALSE)$edge[1:3])
  expect_warning(expect_equal(RandomTree(3, root = "t2"),
                              RandomTree(3, root = 2:3)))
  expect_error(RandomTree(4, root = "not_there"), "No match found for `root`")
  expect_error(RandomTree(4, root = 999), "exceeds number of nodes")
  expect_error(RandomTree(4, root = -1), "must be a positive integer")
    
  expect_warning(expect_equal({
    set.seed(1)
    RandomTree(8, root = NA_integer_)
  }, {
    set.seed(1)
    RandomTree(8, root = FALSE)
  }), "Treating `root = NA` as `FALSE`")
  expect_error(RandomTree(4, nodes = 0), "A tree must contain one .*node")

  expect_warning(RandomTree(4, nodes = 4),
                 "`nodes` higher than number in binary tree")
  expect_warning(RandomTree(4, root = FALSE, nodes = 3),
                 "`nodes` higher than number in binary tree")

  for (nNode in 1:8) {
    expect_equal(RandomTree(10, nodes = nNode)$Nnode, nNode)
  }
  for (nNode in 1:9) {
    expect_equal(RandomTree(10, root = TRUE, nodes = nNode)$Nnode, nNode)
  }
})

test_that("RandomTree(lengths)", {
  expect_equal(RandomTree(0, lengths = 1), ZeroTaxonTree())
  expect_equal(RandomTree(1, lengths = 2), SingleTaxonTree("t1", lengths = 2))
  expect_equal(RandomTree(3, lengths = 2)$edge.length, c(2, 2, 2))
  expect_equal(RandomTree(3, lengths = 1:3, root = FALSE)$edge.length, 1:3)
  expect_equal(RandomTree(3, lengths = 1:3, root = TRUE)$edge.length, c(1:3, 1))
  expect_equal(RandomTree(5, nodes = 2, lengths = 1:3)$edge.length, c(1:3, 1:3))
})

test_that("Small random trees are generated", {
  expect_equal(RandomTree(0, root = FALSE), ZeroTaxonTree())
  expect_equal(RandomTree(0, root = TRUE), ZeroTaxonTree())
  expect_equal(RandomTree(0), ZeroTaxonTree())
  expect_equal(RandomTree("one", root = FALSE), SingleTaxonTree("one"))
  expect_equal(RandomTree("one", root = TRUE), SingleTaxonTree("one"))
  expect_equal(RandomTree(letters[1:2], root = FALSE),
               UnrootTree(PectinateTree(letters[1:2])))
  expect_equal(RandomTree(letters[1:2], root = TRUE),
               UnrootTree(PectinateTree(letters[1:2])))
})

test_that("YuleTree() works", {
  expect_equal(YuleTree(0), ZeroTaxonTree())
  expect_equal(YuleTree("a"), SingleTaxonTree("a"))
  expect_equal(YuleTree(2, addInTurn = TRUE), BalancedTree(2))
  expect_equal(NTip(YuleTree(10)), 10)
  expect_equal(YuleTree(10)$Nnode, 9)
  expect_true(all(TipLabels(YuleTree(10)) %in% TipLabels(10)))
  expect_equal(TipLabels(YuleTree(10, addInTurn = TRUE)), TipLabels(10))
  expect_equal(
    mean(replicate(100, TotalCopheneticIndex(YuleTree(10)))),
    TCIContext(10)$yule.expected,
    tolerance = 0.1
  )
})

test_that("YuleTree() root parameter", {
  expect_equal(
    {set.seed(0); YuleTree(10, root = FALSE)},
    {set.seed(0); UnrootTree(YuleTree(10, root = TRUE))}
  )
  
  expect_warning(expect_equal(
    {set.seed(0); YuleTree(10, root = NA)},
    {set.seed(0); YuleTree(10, root = FALSE)}
  ), "root = NA")
})

test_that("YuleTree(lengths)", {
  expect_equal(YuleTree(0, lengths = 1), ZeroTaxonTree())
  expect_equal(YuleTree(1, lengths = 2), SingleTaxonTree("t1", lengths = 2))
  expect_equal(YuleTree(3, lengths = 2)$edge.length, rep(2, 4))
  expect_equal(YuleTree(3, lengths = 1:3, root = FALSE)$edge.length, 1:3)
  expect_equal(YuleTree(3, lengths = 1:3, root = TRUE)$edge.length, c(1:3, 1))
})

test_that("Hamming() works", {
  dataset <- StringToPhyDat("111100 ???000 ???000 111??? 10??10",
                            letters[1:5], byTaxon = TRUE)
  expected <- c(1/3, 1/3, 0, 1/2,
                0, NaN, 1/2,
                NaN, 1/2,
                1/2)
  expect_equal(as.double(Hamming(dataset, ambig = "NAN")), expected)
  ex <- expected
  ex[is.nan(expected)] <- NA
  expect_equal(as.double(Hamming(dataset, ambig = c("NA", "mean"))), ex)
  ex[is.nan(expected)] <- 0
  expect_equal(as.double(Hamming(dataset, ambig = 0)), ex)
  ex[is.nan(expected)] <- 1
  expect_equal(as.double(Hamming(dataset, ambig = "1")), ex)
  ex[is.nan(expected)] <- mean(expected[!is.nan(expected)])
  expect_equal(as.double(Hamming(dataset, ambig = "mean")), ex)
  ex[is.nan(expected)] <- median(expected[!is.nan(expected)])
  expect_equal(as.double(Hamming(dataset, ambig = "med")), ex)
  expect_error(Hamming(dataset, ambig = "ERROR"))
})

test_that("Hamming() handles inapplicables", {
  dataset <- StringToPhyDat("221100 ---000 ---000 211{-0}?? 10-?10",
                            letters[1:5], byTaxon = TRUE)
  expected <- c(1/3, 1/3, 1/3, 3/4,
                0, NaN, 1/2,
                NaN, 1/2,
                1)
  expect_equal(as.double(Hamming(dataset, ambig = "NaN")), expected)

})

test_that("NJTree() works", {
  a..f <- letters[1:6]
  bal6 <- StringToPhyDat("111100 111000 111000 110000", letters[1:6],
                         byTaxon = FALSE)
  expect_true(all.equal(
    RootTree(NJTree(bal6), a..f[1:3]),
    BalancedTree(letters[c(1:3, 6:4)])
  ))
  expect_equal(NJTree(bal6, edgeLengths = TRUE),
               Preorder(NJTree(bal6, edgeLengths = TRUE)))
  expect_equal(c(0, 1, 2, 1, rep(0, 6)),
               RootTree(NJTree(bal6, edgeLengths = TRUE), 6)$edge.length * 4L)
})

test_that("Constrained NJ trees work", {
  dataset <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 1,
      0, 1, 1, 0, 0, 1), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  constraint <- MatrixToPhyDat(c(a = 0, b = 0, c = 0, d = 0, e = 1, f = 1))
  expect_true(all.equal(read.tree(text = "(a, (d, ((c, b), (e, f))));"),
                        ConstrainedNJ(dataset, constraint)))
  # b == c == f, so these three could be resolved in one of three ways. Drop B.
  expect_true(all.equal(DropTip(NJTree(dataset), "b"),
                        DropTip(ConstrainedNJ(dataset, dataset), "b")))

  expect_true(all.equal(
    KeepTip(ConstrainedNJ(dataset, constraint[3:6]), letters[3:6]),
    BalancedTree(letters[3:6])
  ))
})

test_that("Hamming() fails nicely", {
  expect_error(Hamming(matrix(1:4, 2, 2)))
})
