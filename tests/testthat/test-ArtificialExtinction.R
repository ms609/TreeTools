context("ArtificalExtinction.R")

nTax <- 6L
dataset <- matrix(seq_len(nTax), nrow = nTax, ncol = 6,
                  dimnames = list(LETTERS[1:nTax], 1:6))
dataset[4, 3:5] <- '-'
dataset[5:6, 4:6] <- '?'
dataset[6, 3] <- '?'

test_that("Errors are handled", {

  expect_error(ArtificialExtinction(dataset, 'A', 'E', replaceAmbiguous = 'INVALID'))
  expect_error(ArtificialExtinction(dataset, 'A', 'E', replaceCoded = 'INVALID'))
})

test_that("Replacements ok", {
  expect_equivalent(rep('?', 6),
                    ArtificialExtinction(dataset,
                                         c('A', 'B'), 'E', 'ambig')[1:2, 4:6])

  expect_equivalent(rep('1', 6),
                    ArtificialExtinction(dataset, c('A', 'B'), 'F',
                                         'unif', sampleFrom = 'A')[1:2, 4:6])

  expect_true(all(!'?' == ArtificialExtinction(dataset[-6, ], 1:2, 5, 'freq')[1:2, 4:6]))

  nChar <- 100
  dataset <- rbind(subj = c(rep(2, nChar)),
                   templ = c(rep(3, nChar / 2), rep('?', nChar / 2)))
  suppressWarnings(expect_gt(chisq.test(as.integer(
    ArtEx(dataset, 1, 2, 'binary', 'binary')[1, ]))$p.value, 0.05))
})
