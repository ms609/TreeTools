test_that("Integers fit into integer sizes", {
  INT_MAX <- 2147483647L
  expect_equal(as.integer64(INT_MAX), as.integer64("2147483647"))
  expect_gt(.TT_BASE[1], .TT_BASE[2])
})

test_that("as.phylo.numeric()", {
  expect_true(all.equal(as.phylo(0:2, 6, letters[1:6])[[1]],
                        PectinateTree(letters[c(1, 3:6, 2)])))
  expect_true(all.equal(as.phylo(0, tipLabels = letters[1:6]),
                        PectinateTree(letters[c(1, 3:6, 2)])))
  expect_error(as.phylo(0))
  expect_equal(as.phylo(123, nTip = 8),
               as.phylo(123L, nTip = 8))
  expect_equal(as.phylo(123, nTip = 8),
               as.phylo(as.integer64(123), nTip = 8))
  expect_equal(as.phylo(123:124, nTip = 10),
               as.phylo(as.integer64(123:124), nTip = 10))
  expect_equal(as.phylo(as.integer64(10), tipLabels = paste0("t", 1:8)),
               as.phylo(10, 8))
})

test_that("as.TreeNumber() error handling", {
  expect_error(as.phylo(integer64(1), NULL, NULL))
  expect_error(as.TreeNumber(StarTree(5)))
  expect_error(as.TreeNumber(BalancedTree(52)))  # > 51 leaves not supported
})

test_that("as.TreeNumber()", {
  expect_equal(SingleTaxonTree("t1"), as.phylo(integer64(1), 1, NULL))
  expect_equal(c("Phylogenetic tree number 0 of 105 ",
                 " 6 tips: t1 t2 t3 t4 t5 t6"),
               capture.output(print(as.TreeNumber(as.phylo(105, 6)))))
  expect_equal("Phylogenetic tree number 105 of 6332659870762850625 ",
               capture.output(print(as.TreeNumber(as.phylo(105, 19))))[1])
  # 20+ leaves: no warning, decimal string storage
  expect_no_warning(as.TreeNumber(BalancedTree(20)))
  expect_equal("Phylogenetic tree number 105 of 2.216431e+20 ",
               capture.output(print(as.TreeNumber(as.phylo(105, 20))))[1])
  expect_equal(unlist(lapply(1:3, as.integer64)),
               unlist(as.TreeNumber(as.phylo(1:3, 6))))

  tn6 <- as.TreeNumber(as.phylo(1, 6, letters[1:6]))
  expect_equal(as.phylo(1, 6), as.phylo(tn6, tipLabels = NULL))

  tn16 <- as.TreeNumber(as.phylo(1, 16, letters[1:16]))
  expect_equal(as.phylo(1, 16), as.phylo(tn16, tipLabels = NULL))
  expect_equal(as.TreeNumber(as.MixedBase(tn16)), as.TreeNumber(tn16))

  expect_error(as.phylo(as.integer64(2) ^ 62, 16), "too large ")
  bigNumber <- as.integer64(2) ^ 61 + 1
  # Not necessarily equal to as.integer64(2 ^ 61 + 1), due to numeric conversion
  tn16 <- as.TreeNumber(as.phylo(bigNumber, 16, letters[1:16]))
  expect_equal(as.phylo(bigNumber, 16), as.phylo(tn16, tipLabels = NULL))
  # Rounding:
  expect_identical(as.integer64(as.TreeNumber(as.phylo(bigNumber, 16))) - 1,
                   as.integer64(as.TreeNumber(as.phylo(bigNumber - 1, 16))))
               
  expect_equal(as.TreeNumber(as.MixedBase(tn16)), tn16)

  # Large trees: character-backed, round-trips correctly.
  # as.phylo.TreeNumber returns Preorder(RootTree(x, 1)), so compare against that.
  for (nTip in c(20L, 35L, 51L)) {
    tree <- BalancedTree(nTip)
    tn <- as.TreeNumber(tree)
    expect_s3_class(tn, "TreeNumber")
    expect_true(inherits(tn, "character"))
    expect_false(inherits(tn, "integer64"))
    expect_equal(Preorder(RootTree(tree, 1)), as.phylo(tn, tipLabels = NULL))
  }
})

test_that("as.MixedBase()", {
  expect_error(as.MixedBase(RandomTree(100000)),
               "Too many \\w+ for mixed base")
  expect_equal(as.integer(as.TreeNumber(as.phylo(16, 16))), 16)
  nTip <- 9
  expect_equal(as.integer(as.MixedBase(as.TreeNumber(as.phylo(16, nTip)))),
               c(rep(0, nTip - 3 - 3), 1, 0, 1))
  
  expect_equal(as.MixedBase(as.phylo(0, 6)),
               structure(integer(3),
                         nTip = 6,
                         tip.label = paste0("t", 1:6),
                         binary = TRUE,
                         class = "MixedBase"))
  
  expect_equal(as.MixedBase(as.phylo(105, 20)),
               structure(tabulate(17 - 3, 17),
                         nTip = 20,
                         tip.label = paste0("t", 1:20),
                         binary = TRUE,
                         class = "MixedBase"))
  
  expect_equal(unlist(lapply(1:3, function(n) as.MixedBase(as.phylo(n, 6)))),
               unlist(as.MixedBase(as.phylo(1:3, 6))))
  
  mb3 <- as.MixedBase(as.phylo(1, 3, letters[1:3]))
  expect_equal(as.phylo(1, 3), as.phylo(mb3, tipLabels = NULL))
  
  mb6 <- as.MixedBase(as.phylo(1, 6, letters[1:6]))
  expect_equal(as.phylo(1, 6), as.phylo(mb6, tipLabels = NULL))
  
  mb16 <- as.TreeNumber(as.phylo(1337, 16, letters[1:16]))
  expect_equal(as.phylo(1337, 16), as.phylo(mb16, tipLabels = NULL))
  
  expect_error(as.MixedBase(42, TipLabels(2)))
  expect_true(as.MixedBase(as.phylo(1337, 16)) == as.MixedBase(1337, 16))
  expect_true(as.MixedBase(44) > as.MixedBase(as.MixedBase(42)))
  expect_true(as.MixedBase(44, 10) < as.MixedBase(42, 11))
  bigNo <- as.integer64(2) ^ 62 + 1337
  expect_equal(sum(as.integer(as.MixedBase(bigNo)) * .TT_BASE), bigNo)
})

test_that("as.MixedBase() supports larger trees", {
  expect_equal(tail(as.MixedBase(BalancedTree(100))),
               c(6, 10, 8, 3, 4, 1) # From observation, not calculation
  )
})

test_that("is.TreeNumber()", {
  expect_equal(is.TreeNumber(as.integer64(5)), FALSE)
  expect_equal(is.TreeNumber(as.TreeNumber(BalancedTree(5))), TRUE)
})

test_that("tree_number.h add_small carry propagation", {
  # The tree number 2^64 + 1 = "18446744073709551617" has packed chunks [4, 8, 5].
  # Reconstructing it via packed_to_tree_num() calls add_small(5) after mul_small()
  # leaves w[0] = 2^64 - 4, so w[0] + 5 overflows, triggering the carry-propagation
  # loop (lines 86-88 of inst/include/TreeTools/tree_number.h).
  # NUnrooted(20) > 2^64, so this is a valid 20-leaf tree number.
  tn <- as.TreeNumber("18446744073709551617", nTip = 20L)
  expect_true(inherits(tn, "character"))
  expect_equal(tn, as.TreeNumber(as.phylo(tn, tipLabels = NULL)))
})

test_that(".decimal_to_chunks() handles zero and empty input", {
  # Calling as.phylo() on a character-backed TreeNumber of "0" routes through
  # .decimal_to_chunks("0"), hitting the early-return branch.
  tn_zero <- as.TreeNumber("0", nTip = 20L)
  expect_equal(as.phylo(tn_zero, tipLabels = NULL),
               as.phylo(0L, nTip = 20L))

  # Direct internal check: both "0" and "" return 0L immediately.
  expect_equal(.decimal_to_chunks("0"), 0L)
  expect_equal(.decimal_to_chunks(""), 0L)

  # The else-branch is always taken when input is non-zero;
  # the `0L` branch is unreachable (while loop always appends >= 1 element).
  expect_equal(.decimal_to_chunks("1"), 1L)
})

test_that("as.TreeNumber.character() creates character-backed result for nTip > .TT_MAX_TIP", {
  tn <- as.TreeNumber("12345678901234", nTip = 25L)
  expect_s3_class(tn, "TreeNumber")
  expect_true(inherits(tn, "character"))
  expect_false(inherits(tn, "integer64"))
  expect_equal(attr(tn, "nTip"), 25L)
})

test_that("as.TreeNumber.MixedBase() errors for trees with > .TT_MAX_TIP leaves", {
  mb <- as.MixedBase(BalancedTree(20))
  expect_error(as.TreeNumber(mb), "MixedBase.*TreeNumber")
  tn <- as.TreeNumber(BalancedTree(20))
  expect_error(as.MixedBase(tn), "TreeNumber.*MixedBase")
})

test_that("as.integer64.TreeNumber() for character-backed TreeNumbers", {
  # converting a large (nTip > 19) character-backed TreeNumber
  # to integer64 is impossible; should error.
  tn_large <- as.TreeNumber(BalancedTree(20))
  expect_error(as.integer64(tn_large), "integer64")

  # a character-backed TreeNumber with nTip <= 19 (atypical but valid
  # as a manually-constructed object) can be converted exactly to integer64.
  tn_small_char <- structure("10", nTip = 6L, tip.label = paste0("t", 1:6),
                             class = c("TreeNumber", "character"))
  expect_equal(as.integer64(tn_small_char), as.integer64(10))
})
