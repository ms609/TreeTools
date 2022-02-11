test_that("Node supports calculated correctly", {
  treeSample <- list(
    correct = ape::read.tree(text = "((((((A,B),C),D),E),F),out);"),
    swapFE  = ape::read.tree(text = "((((((A,B),C),D),F),E),out);"),
    DEClade = ape::read.tree(text = "(((((A,B),C),(D,E)),F),out);"),
    swapBC  = ape::read.tree(text = "((((((A,C),B),D),E),F),out);"),
    DbyA    = ape::read.tree(text = "((((((A,D),C),B),E),F,G),out);")
  )
  expect_equal(c('10' = 4, '11' = 4, '12' = 4, '13' = 3),
               SplitFrequency(treeSample$correct, treeSample))

  # Internal nodes on each side of root
  balanced <- ape::read.tree(text="((D, (E, (F, out))), (C, (A, B)));")
  freq <- SplitFrequency(balanced, treeSample)
  expect_equal(freq,
               c("9" = 4, "10" = 4, "11" = 4, "12" = 4, "13" = 3)[names(freq)])

})

test_that("Node support colours consistent", {
  expect_equal('red', SupportColour(NA))
  expect_equal(c('#ffffff00', 'red'), SupportColour(1:2, show1 = FALSE))
  expect_equal('red', SupportColor(-2)) # Check alternative spelling
  expect_equal('#ffffff00', SupportColour(1, show1 = FALSE))
  expect_equal(c('oor', '1', '34', '67', '101'),
               SupportColour((-1):3 / 3, scale = 1:101, outOfRange = 'oor'))
})

test_that("SplitFrequency() handles four-split trees", {
  trees <- AddTipEverywhere(BalancedTree(3))
  trees <- c(trees[1], trees)
  expect_equal(c('7' = 2L), SplitFrequency(trees[[1]], trees))
})

test_that("LabelSplits()", {
  expect_error(LabelSplits(BalancedTree(8), 1:8))
  skip_if_not_installed('vdiffr', minimum_version = "1.0.0")
  skip_if(packageVersion("graphics") < "4.1.0")
  vdiffr::expect_doppelganger('LabelSplits()',  function() {
    tree <- BalancedTree(9)
    plot(tree)
    labs <- letters[6:1]
    names(labs) <- rev(names(as.Splits(tree)))
    LabelSplits(tree, labs, frame = 'circ', cex = 2, bg = 'orange')
  })
  vdiffr::expect_doppelganger('LabelSplits()-nameless', function() {
    tree <- BalancedTree(9)
    plot(tree)
    LabelSplits(tree, bg = 'orange')
    expect_warning(LabelSplits(BalancedTree(9), setNames(letters[11:16], 1:6)))
  })
  vdiffr::expect_doppelganger('LabelSplits()-names', function() {
    tree <- BalancedTree(9)
    plot(tree)
    labs <- letters[1:6]
    LabelSplits(tree, labs, bg = 'orange')
  })
})
