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
  expectation <- matrix(rep('?', 6), 2, 3, dimnames = list(c('A', 'B'), 4:6))
  expect_equal(expectation,
               ArtificialExtinction(dataset, subject = c('A', 'B'),
                                    template = 'E',
                                    replaceAmbiguous = 'ambig')[1:2, 4:6])

  expectation[] <- 1
  expect_equal(expectation,
               ArtificialExtinction(dataset, c('A', 'B'), 'F',
                                    'unif', sampleFrom = 'A')[1:2, 4:6])

  expect_true(all(!'?' == ArtificialExtinction(dataset[-6, ], 1:2, 5, 'freq')[1:2, 4:6]))


  expect_true(all(ArtificialExtinction(dataset, subject = 'E', template = 'F',
                                       replaceAmbiguous = 'binary')['E', 3:6] %in% 0:1))

  expect_equal(setNames(rep('?', 3), 4:6),
               ArtificialExtinction(dataset, subject = 'E', template = 'F',
                                    replaceAmbiguous = 'binary',
                                    replaceAll = FALSE)['E', 4:6])

  nChar <- 100
  dataset <- rbind(subj = c(rep(2, nChar)),
                   templ = c(rep(3, nChar / 2), rep('?', nChar / 2)))
  suppressWarnings(expect_gt(chisq.test(as.integer(
    ArtEx(dataset, 1, 2, 'binary', 'binary')[1, ]))$p.value, 0.05))
})
