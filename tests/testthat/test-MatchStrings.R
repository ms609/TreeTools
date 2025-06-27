test_that("MatchStrings()", {
  list1 <- "Mouse"
  list2 <-c("cat", "mouse", "dog")
  expect_error(MatchStrings("maus", list1, stop),
               "Could not find 'maus' in list1.  Did you mean 'Mouse'?")
  expect_warning(
    expect_equal(
      MatchStrings(c("cat", "mouse", "dawg", "maus", "non-existent taxon"),
                   list2, warning),
      c("cat", "mouse")),
    "Could not find 'dawg', 'maus', 'no.+on' in list2.+ mean 'dog', 'mouse'\\?"
  )
})
