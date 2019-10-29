context("Support.R")

test_that("Node supports calculated correctly", {
  treeSample <- list(
    correct = ape::read.tree(text = "((((((A,B),C),D),E),F),out);"),
    swapFE  = ape::read.tree(text = "((((((A,B),C),D),F),E),out);"),
    DEClade = ape::read.tree(text = "(((((A,B),C),(D,E)),F),out);"),
    swapBC  = ape::read.tree(text = "((((((A,C),B),D),E),F),out);"),
    DbyA    = ape::read.tree(text = "((((((A,D),C),B),E),F),out);")
  )
  expect_equal(c('n10'=4, 'n11'=4, 'n12'=4, 'n13'=3),
               SplitFrequency(treeSample$correct, treeSample))

  # Internal nodes on each side of root
  balanced <- ape::read.tree(text="((D, (E, (F, out))), (C, (A, B)));")
  expect_equal(c('n10'=4, 'n11'=4, 'n12'=4, 'n13'=3),
               SplitFrequency(balanced, treeSample))

})

test_that("Node support colours consistent", {
  expect_equal('red', SupportColour(NA))
  expect_equal('red', SupportColour(2))
  expect_equal('red', SupportColor(-2)) # Check alternative spelling
  expect_equal('#ffffff00', SupportColour(1, show1=FALSE))
})
