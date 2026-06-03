test_that("Consensus() errors", {
  expect_error(Consensus(trees), " a list")
  
  bal8 <- BalancedTree(8)
  oneLeaf <- Consensus(list(DropTip(bal8, 1:7))[c(1, 1, 1)])
  expect_equal(class(oneLeaf), "phylo")
  expect_equal(oneLeaf$Nnode, 0)
  expect_identical(Consensus(list(DropTip(bal8, 1:8))[c(1, 1, 1)]),
                   DropTip(bal8, 1:8))
  expect_warning(expect_true(all.equal(
    Consensus(list(PectinateTree(6), PectinateTree(8))),
    PectinateTree(6)
  )), regexp = "Tree sizes")

  halfTree <- CollapseNode(bal8, 10:12)
  expect_equal(Consensus(halfTree), halfTree)
  expect_equal(Consensus(c(halfTree)), halfTree)
  expect_equal(Consensus(list(halfTree)), halfTree)
  
  # Test that trees larger than heap limit fail gracefully
  expect_error(
    Consensus(c(BalancedTree(100001), PectinateTree(100001))),
    "too many leaves.*100000"
  )
  
  largeTree <- BalancedTree(33333)
  consensus_large <- Consensus(c(largeTree, largeTree))
  expect_equal(NTip(consensus_large), 33333)
})

test_that("Consensus()", {

  ApeTest <- function(tr, p = 1) {
    if (!expect_true(all.equal(RootTree(Consensus(tr, p = p), 1),
                               RootTree(ape::consensus(tr, p = p), 1)))) {
      dput(RootTree(Consensus(tr, p = p), 1)$edge)
      dput(RootTree(ape::consensus(tr, p = p), 1)$edge)
    }
  }

  trees <- list(BalancedTree(8), PectinateTree(8))[c(1, 1, 1, 1, 2, 2, 2)]
  trees <- RenumberTips(trees, trees[[1]])
  expect_error(Consensus(trees, 2), regexp = "`p`") # p too high
  expect_error(Consensus(trees, 0.2), regexp = "`p`") # p too low

  ApeTest(trees)
  ApeTest(trees, 0.5)
  expect_equal(Consensus(trees, 0.999), Consensus(trees, 1))

  # Even-n exact tie: two conflicting splits each in exactly 50% of trees.
  # Strict majority keeps a split only when count > n / 2 (thresh =
  # floor(n / 2) + 1), so both conflicting splits are dropped and the consensus
  # is unresolved there -- always a valid tree.  The "exact == hashed" test
  # below exercises this fixture too, but only against the other internal
  # counter; pin it to an independent oracle (ape) here.
  tie <- c(rep(list(BalancedTree(8)), 2L), rep(list(PectinateTree(8)), 2L))
  tie <- RenumberTips(tie, tie[[1]])
  ApeTest(tie, 0.5)

  # p > 0.5 boundary: a split in exactly p * n trees is retained ("p or more"),
  # matching ape.  Clade {t1, t2} occurs in 2 of 3 trees (= 2/3 exactly).
  thirds <- lapply(c("(((t1, t2), t3), t4, t5, t6);",
                     "(((t1, t2), t3), t4, t5, t6);",
                     "((t1, (t2, t3)), t4, t5, t6);"),
                   function(x) ape::read.tree(text = x))
  ApeTest(thirds, 2 / 3)

  ApeTest(as.phylo(0:2, 8))
  ApeTest(as.phylo(0:250, 8))
  ApeTest(as.phylo(0:250, 80))
  
  trees <- list(ape::read.tree(text = "((a, b), (c, d));"),
                ape::read.tree(text = "((a, c), (b, d));"))
  expect_equal(Consensus(trees), Preorder(StarTree(letters[1:4])))
})

test_that("Consensus() ignores edge lengths, node labels, and TipLabel", {
  base <- list(BalancedTree(8), PectinateTree(8))[c(1, 1, 1, 1, 2, 2, 2)]
  base <- structure(RenumberTips(base, base[[1]]), class = "multiPhylo")
  withEL <- structure(lapply(base, function(t) {
    t[["edge.length"]] <- seq_len(nrow(t[["edge"]])); t
  }), class = "multiPhylo")
  withNL <- structure(lapply(base, function(t) {
    t[["node.label"]] <- paste0("n", seq_len(t[["Nnode"]])); t
  }), class = "multiPhylo")
  labelled <- ape::.compressTipLabel(base)
  expect_true(is.null(.subset2(labelled, 1)[["tip.label"]]))  # genuinely labelled

  for (p in c(0.5, 2 / 3, 1)) for (check in c(TRUE, FALSE)) {
    ref <- Consensus(base, p, check.labels = check)
    expect_equal(Consensus(withEL, p, check.labels = check), ref)
    expect_equal(Consensus(withNL, p, check.labels = check), ref)
    expect_equal(Consensus(labelled, p, check.labels = check), ref)
  }
  # Result carries no inherited metadata.
  out <- Consensus(withEL, 0.5)
  expect_null(out[["edge.length"]])
  expect_null(out[["node.label"]])

  # Differing tree sizes exercise the KeepTip repeat-loop; metadata must not
  # change its (warned) result either.
  mixed <- list(BalancedTree(8), PectinateTree(8), DropTip(BalancedTree(8), 1:2))
  mixedEL <- lapply(mixed, function(t) {
    t[["edge.length"]] <- seq_len(nrow(t[["edge"]])); t
  })
  expect_equal(suppressWarnings(Consensus(mixedEL, 0.5)),
               suppressWarnings(Consensus(mixed, 0.5)))
})

test_that("Consensus() handles large sets of trees", {
  oneTree <- as.phylo(0, 13)
  expect_equal(
    Consensus(lapply(1:33000, function(x) oneTree)),
    oneTree
  )
  
  manyTrees <- lapply(1:33000, as.phylo, 13) # More than int16_max
  expect_equal(
    Consensus(c(manyTrees[1:16500], lapply(1:16500, function(x) oneTree)),
              p = 0.5),
    oneTree
  )
  
  skip_on_cran() # Slow!
  skip_if(!isTRUE(options("runSlowTests")))
  
  expect_true(all.equal(
    Consensus(manyTrees),
    read.tree(text = "((((((t2,t13),t12),t11),t10),t3,t4,t5,t6,t7,t8,t9),t1);")
    # i.e. write.tree(ape::consensus(manyTrees))
  ))
})

test_that("Consensus() exact and hashed counts agree", {
  # The hashed (default) and exact (opt-in) split counts must yield identical
  # consensus trees; this also guards the shared counting core.
  skip_if_not_installed("ape")
  set.seed(1)
  forests <- list(
    balPec = list(BalancedTree(8), PectinateTree(8))[c(1, 1, 1, 1, 2, 2, 2)],
    starlike = list(ape::read.tree(text = "((a, b), (c, d));"),
                    ape::read.tree(text = "((a, c), (b, d));")),
    tie = c(rep(list(BalancedTree(8)), 2L), rep(list(PectinateTree(8)), 2L)),
    rand12 = lapply(1:7, function(i) ape::rtree(12, br = NULL)),
    rand9  = lapply(1:20, function(i) ape::rtree(9, br = NULL))
  )
  for (f in forests) {
    for (p in c(0.5, 2 / 3, 1)) {
      hashed <- Consensus(f, p = p)
      exact  <- Consensus(f, p = p, hash = FALSE)
      expect_true(isTRUE(all.equal(RootTree(hashed, 1), RootTree(exact, 1))))
    }
  }
})

test_that("Consensus() handles non-preorder trees", {
  trees <- ape::read.nexus(test_path("testdata", "nonPreCons.nex"))
  expect_equal(Consensus(trees)$Nnode, 3)
  expect_equal(Consensus(trees), Preorder(trees[[1]]))

  # Consensus() now preorders only trees[[1]] (consensus_tree() preorders each
  # tree internally, #168 fixed 2026-06). Feeding genuinely edge-shuffled trees
  # must give the same result as preordering the whole forest first.
  shuffle <- function(t) {
    t2 <- RootTree(t, sample(NTip(t), 1L))
    e  <- t2$edge
    t2$edge <- e[sample(nrow(e)), , drop = FALSE]
    attr(t2, "order") <- NULL
    t2
  }
  set.seed(3L)
  base <- ape::rtree(10L, br = NULL)
  forest <- RenumberTips(
    lapply(1:6, function(i) ape::rtree(10L, br = NULL, tip.label = base$tip.label)),
    base)
  nonpre <- lapply(forest, shuffle)
  for (p in c(0.5, 2 / 3, 1)) {
    expect_equal(Consensus(nonpre, p = p, check.labels = FALSE),
                 Consensus(Preorder(nonpre), p = p, check.labels = FALSE))
  }
})

test_that("consensus_tree() is robust to non-preorder input", {
  # Regression for #168 / root_on_node bug: ClusterTable used to return
  # wrong splits or segfault when input edges were not in preorder.
  # Non-preorder is now handled internally; results must equal the preordered run.

  # fixture: cladewise (non-preorder), differing roots
  raw <- lapply(ape::read.nexus(test_path("testdata", "nonPreCons.nex")),
                identity)  # plain list so each phylo carries its own tip.label
  tr4 <- lapply(raw, function(t) RenumberTips(t, raw[[1]]))
  pre4 <- Preorder(tr4)
  expect_equal(TreeTools:::consensus_tree(tr4,  1, FALSE),
               TreeTools:::consensus_tree(pre4, 1, FALSE))
  expect_equal(nrow(TreeTools:::consensus_tree(tr4, 1, FALSE)), 1L)

  # edge-shuffled trees: random root + permuted edge rows
  set.seed(1L)
  base <- ape::rtree(8L, br = NULL)
  shuffle <- function(t) {
    t2 <- RootTree(t, sample(NTip(t), 1L))
    e  <- t2$edge
    t2$edge <- e[sample(nrow(e)), , drop = FALSE]
    attr(t2, "order") <- NULL
    t2
  }
  tr8 <- RenumberTips(lapply(seq_len(4L), function(i) shuffle(base)), base)
  pre8 <- Preorder(tr8)

  expect_equal(TreeTools:::consensus_tree(tr8, 1,   FALSE),
               TreeTools:::consensus_tree(pre8, 1,   FALSE))
  expect_equal(TreeTools:::consensus_tree(tr8, 0.5, TRUE),
               TreeTools:::consensus_tree(pre8, 0.5, TRUE))
  expect_equal(TreeTools:::consensus_tree(tr8, 0.5, FALSE),
               TreeTools:::consensus_tree(pre8, 0.5, FALSE))
})

test_that("split_frequencies() is robust to non-preorder input", {
  # Same ClusterTable hazard as consensus_tree() — both paths (hashed/exact).
  set.seed(2L)
  base <- ape::rtree(8L, br = NULL)
  shuffle <- function(t) {
    t2 <- RootTree(t, sample(NTip(t), 1L))
    e  <- t2$edge
    t2$edge <- e[sample(nrow(e)), , drop = FALSE]
    attr(t2, "order") <- NULL
    t2
  }
  tr <- RenumberTips(lapply(seq_len(4L), function(i) shuffle(base)), base)
  pre <- Preorder(tr)

  sf_non  <- TreeTools:::split_frequencies(tr,  exact = FALSE)
  sf_pre  <- TreeTools:::split_frequencies(pre, exact = FALSE)
  sfx_non <- TreeTools:::split_frequencies(tr,  exact = TRUE)
  sfx_pre <- TreeTools:::split_frequencies(pre, exact = TRUE)

  expect_equal(sort(sf_non$counts),  sort(sf_pre$counts))
  expect_equal(sort(sfx_non$counts), sort(sfx_pre$counts))
})

test_that("ConsensusWithout() is robust", {
  tr <- as.phylo(0:6, 16)
  expect_identical(ConsensusWithout(tr, paste0("t", 1:16)),
                   DropTip(tr[[1]], 1:16))
  expect_identical(ConsensusWithout(tr, paste0("t", 1:4)),
                   Consensus(DropTip(tr, 1:4)))

  expect_true(all.equal(BalancedTree(8), ConsensusWithout(BalancedTree(8))))
  expect_true(all.equal(BalancedTree(4),
                        ConsensusWithout(BalancedTree(8), paste0("t", 5:8))))
  balAndPec <- list(BalancedTree(8), PectinateTree(8))
  t25 <- paste0("t", c(2:5))
  expect_true(all.equal(PectinateTree(paste0("t", c("1", 6:8))),
                        ConsensusWithout(balAndPec, t25)))
  expect_true(all.equal(
    ConsensusWithout(structure(balAndPec, class = "multiPhylo"), t25),
    ConsensusWithout(balAndPec, t25))
  )

  nasty <- structure(list(edge = structure(
    c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
      5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
    .Dim = c(12, 2)),
    Nnode = 5L,
    tip.label = letters[1:8]),
    class = "phylo") # Danger: Do not plot!
  expect_true(all.equal(Preorder(nasty), ConsensusWithout(nasty)))
  expect_true(all.equal(DropTip(nasty, 2), ConsensusWithout(nasty, "b")))

})
