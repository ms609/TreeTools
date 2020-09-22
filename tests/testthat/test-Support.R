context("Support.R")

test_that("Node supports calculated correctly", {
  treeSample <- list(
    correct = ape::read.tree(text = "((((((A,B),C),D),E),F),out);"),
    swapFE  = ape::read.tree(text = "((((((A,B),C),D),F),E),out);"),
    DEClade = ape::read.tree(text = "(((((A,B),C),(D,E)),F),out);"),
    swapBC  = ape::read.tree(text = "((((((A,C),B),D),E),F),out);"),
    DbyA    = ape::read.tree(text = "((((((A,D),C),B),E),F),out);")
  )
  expect_equal(c('10'=4, '11'=4, '12'=4, '13'=3),
               SplitFrequency(treeSample$correct, treeSample))

  # Internal nodes on each side of root
  balanced <- ape::read.tree(text="((D, (E, (F, out))), (C, (A, B)));")
  expect_equal(c('10'=4, '11'=4, '12'=4, '13'=3),
               SplitFrequency(balanced, treeSample))

})

test_that("Node support colours consistent", {
  expect_equal('red', SupportColour(NA))
  expect_equal(c('#ffffff00', 'red'), SupportColour(1:2, show1 = FALSE))
  expect_equal('red', SupportColor(-2)) # Check alternative spelling
  expect_equal('#ffffff00', SupportColour(1, show1 = FALSE))
  expect_equal(c('oor', '1', '34', '67', '101'),
               SupportColour((-1):3 / 3, scale = 1:101, outOfRange = 'oor'))
})
