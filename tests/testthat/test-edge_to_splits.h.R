Hash <- function(tree) cpp_topology_hash(tree[["edge"]], NTip(tree))

test_that("topology_hash() ignores root position, edge order and node labels", {
  for (nTip in c(6, 70)) {
    tree <- RandomTree(nTip, root = FALSE)
    rooted <- c(
      lapply(c(1, 2, nTip), function(tip) RootTree(tree, tip)),
      list(RootTree(tree, c(2, 3)), UnrootTree(RootTree(tree, 5)),
           Cladewise(RootTree(tree, 4)))
    )
    expect_equal(vapply(rooted, Hash, double(1)),
                 rep(Hash(tree), length(rooted)))

    # Relabel internal nodes without moving any child ahead of its parent
    relabelled <- tree
    internal <- NTip(tree) + seq_len(tree[["Nnode"]])
    newLabel <- c(seq_len(NTip(tree)), NTip(tree) + 1,
                  NTip(tree) + 1 + rev(seq_len(tree[["Nnode"]] - 1)))
    relabelled[["edge"]][] <- newLabel[tree[["edge"]]]
    expect_equal(Hash(relabelled), Hash(tree))
  }
})

test_that("topology_hash() distinguishes topologies", {
  trees <- as.phylo(0:104, 6)
  expect_equal(length(unique(vapply(trees, Hash, double(1)))), 105)

  big <- BalancedTree(70)
  swapped <- big
  swapped[["edge"]][match(c(1, 70), big[["edge"]][, 2]), 2] <- c(70, 1)
  expect_false(Hash(big) == Hash(swapped))
  expect_false(Hash(big) == Hash(CollapseNode(big, 80)))
})

test_that("topology_hash() handles trees without informative splits", {
  empty <- Hash(StarTree(8))
  expect_equal(Hash(BalancedTree(3)), empty)
  expect_equal(Hash(RootTree(BalancedTree(3), 1)), empty)
})
