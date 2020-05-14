context('phylo.R')

test_that('AddTipEverywhere correct', {
  backbone <- PectinateTree(5)

  expect_equal(7L, length(AddTipEverywhere(backbone, includeRoot = FALSE)))
  expect_equal(7L, length(AddTipEverywhere(unroot(backbone), includeRoot = FALSE)))
  expect_equal(9L, length(AddTipEverywhere(backbone, includeRoot = TRUE)))
  expect_equal(5L, length(AddTipEverywhere(CollapseNode(backbone, 7:9))))

})

test_that('AddTipEverywhere handles nasty tree', {
  nasty <- structure(list(edge = structure(
    c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
      5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
    .Dim = c(12, 2)),
    Nnode = 5L,
    tip.label = letters[1:8]),
    class = 'phylo') # Danger: Do not plot!
  expect_equal(AddTipEverywhere(Preorder(nasty)),
               lapply(AddTipEverywhere(nasty), Preorder))
})
