# Tests for the internal fastmatch shim (R/fastmatch.R, .onLoad in R/zzz.R).
#
# fastmatch is Suggested, not Imported (it cannot compile to WebAssembly), so
# TreeTools must fall back to base match()/%in% when it is absent. These tests
# guarantee that .FastMatch() and %fin% return results identical to fastmatch for
# every relevant input, on both code paths.

test_that(".FastMatch() and %fin% equal base match across edge cases", {
  # The base fallback must return results identical to fastmatch for every
  # input type that TreeTools feeds to these functions: integer / double node
  # vectors (e.g. colSums() output) and character tip labels, with NA /
  # no-match / empty / duplicate inputs. Asserted against base match() directly,
  # so it must hold whether or not fastmatch is the active backend.
  #
  # Factors and logicals are deliberately excluded, as TreeTools never matches
  # them via these helpers and fastmatch treats them specially: fmatch() hashes
  # a factor's integer codes (whereas base matches by label, so they genuinely
  # disagree -- fmatch returns NA), and warns + falls back to match() for
  # logicals. Asserting identity on those would be a false equivalence / noise.
  fin <- function(x, table) match(x, table, nomatch = 0L) > 0L

  cases <- list(
    integer    = list(x = c(3L, 1L, 9L, 1L), table = c(1L, 2L, 3L, 4L)),
    no_match   = list(x = c(5L, 6L),         table = c(1L, 2L, 3L)),
    with_na    = list(x = c(2L, NA, 1L),     table = c(NA, 1L, 2L)),
    empty_x    = list(x = integer(0),        table = c(1L, 2L)),
    empty_tab  = list(x = c(1L, 2L),         table = integer(0)),
    double     = list(x = c(2.5, 1.0),       table = c(1.0, 2.5, 3.0)),
    character  = list(x = c("b", "z", "a"),  table = letters[1:5]),
    char_na    = list(x = c("a", NA, "q"),   table = c("a", "b", NA)),
    duplicates = list(x = c(1L, 1L, 2L),     table = c(2L, 1L, 1L))
  )

  for (nm in names(cases)) {
    x <- cases[[nm]][["x"]]
    table <- cases[[nm]][["table"]]
    expect_identical(.FastMatch(x, table), match(x, table), info = nm)
    expect_identical(x %fin% table, fin(x, table), info = nm)
  }

  # nomatch argument is honoured identically
  expect_identical(.FastMatch(c(5L, 1L), 1:3, nomatch = -1L),
                   match(c(5L, 1L), 1:3, nomatch = -1L))
})

test_that(".FastMatch and %fin% resolve to fastmatch when available", {
  if (requireNamespace("fastmatch", quietly = TRUE) &&
      isTRUE(getOption("TreeTools.fastmatch", TRUE))) {
    # .onLoad rebinds the internal names to fastmatch's own functions, so they
    # are the *same object* -- proving zero wrapper/dispatch overhead.
    expect_identical(TreeTools:::.FastMatch, fastmatch::fmatch)
    expect_identical(get("%fin%", envir = asNamespace("TreeTools")),
                     fastmatch::`%fin%`)
  } else {
    expect_identical(TreeTools:::.FastMatch, base::match)
  }
})

test_that("Functions using %fin% give correct results", {
  # Exercise representative call sites end-to-end (results must not depend on
  # which match backend is active).
  bal <- BalancedTree(8)
  expect_equal(NTip(bal[["edge"]]), 8L)
  expect_equal(RootNode(bal[["edge"]]), 9L)
  expect_equal(NodeNumbers(bal), 9:15)

  collapsed <- CollapseNode(bal, 11)
  expect_true(inherits(collapsed, "phylo"))
  expect_equal(NTip(collapsed), 8L)
})
