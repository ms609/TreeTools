context("ArtificalExtinction.R")

nTax <- 6L
dataset <- matrix(seq_len(nTax), nrow = nTax, ncol = 6,
                  dimnames = list(LETTERS[1:nTax], 1:6))
dataset[4, 3:5] <- '-'
dataset[5:6, 4:6] <- '?'
dataset[6, 3] <- '?'

test_that("Errors are handled", {

  expect_error(ArtificialExtinction(dataset, 'E', 'A', replacement = 'INVALID'))
})

test_that("Replacements ok", {
  expect_equivalent(rep('?', 6),
                    ArtificialExtinction(dataset,
                                         'E', c('A', 'B'), '?')[1:2, 4:6])

  expect_equivalent(rep('1', 6),
                    ArtificialExtinction(dataset, 'F', c('A', 'B'),
                                         'unif', sampleFrom = 'A')[1:2, 4:6])

  expect_true(all(!'?' == ArtificialExtinction(dataset[-6, ], 5, 1:2, 'freq')[1:2, 4:6]))

})
