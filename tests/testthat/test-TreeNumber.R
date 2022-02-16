test_that("Integers fit into integer sizes", {
  INT_MAX <- 2147483647L
  expect_equal(as.integer64(INT_MAX), as.integer64("2147483647"))
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
  expect_equal(as.phylo(as.integer64(10), tipLabels = paste0('t', 1:8)),
               as.phylo(10, 8))
})

test_that("as.TreeNumber() error handling", {
  expect_error(as.phylo(integer64(1), NULL, NULL))
  expect_error(as.TreeNumber(StarTree(5)))
  expect_warning(as.TreeNumber(BalancedTree(20)))
})

test_that("as.TreeNumber()", {
  expect_equal(SingleTaxonTree('t1'), as.phylo(integer64(1), 1, NULL))
  expect_equal(c("Phylogenetic tree number 0 of 105 ",
                 " 6 tips: t1 t2 t3 t4 t5 t6"),
               capture.output(print(as.TreeNumber(as.phylo(105, 6)))))
  expect_equal("Phylogenetic tree number 105 of 6332659870762850625 ",
               capture.output(print(as.TreeNumber(as.phylo(105, 19))))[1])
  expect_warning(
    expect_equal("Phylogenetic tree number 105 of 2.216431e+20 ",
                 capture.output(print(as.TreeNumber(as.phylo(105, 20))))[1])
    )
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
  
})

test_that("as.MixedBase()", {
  expect_equal(as.MixedBase(as.phylo(0, 6)),
               structure(integer(3),
                         nTip = 6,
                         tip.label = paste0('t', 1:6),
                         binary = TRUE,
                         class = "MixedBase"))
  
  expect_equal(as.MixedBase(as.phylo(105, 20)),
               structure(tabulate(17 - 3, 17),
                         nTip = 20,
                         tip.label = paste0('t', 1:20),
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
  expect_equal(sum(as.integer(as.MixedBase(bigNo)) * rev(.TT_BASE)), bigNo)
  
})
