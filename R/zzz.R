.onLoad <- function(libname, pkgname) {
  # Resolve the fast-matching backend once, avoiding per-call dispatch overhead.
  # fastmatch is Suggested rather than Imported because it cannot compile to
  # WebAssembly (see R/fastmatch.R). When it is available -- and not disabled via
  # options(TreeTools.fastmatch = FALSE), used to exercise the fallback path in
  # tests -- rebind the internal `.FastMatch` and `%fin%` to fastmatch's own
  # functions, so calls reach them with no wrapper. Otherwise the base
  # equivalents from R/fastmatch.R remain in force.
  if (isTRUE(getOption("TreeTools.fastmatch", TRUE)) &&
      requireNamespace("fastmatch", quietly = TRUE)) {
    ns <- asNamespace(pkgname)
    assign(".FastMatch", fastmatch::fmatch, envir = ns)
    assign("%fin%", fastmatch::`%fin%`, envir = ns)
  }
}

.onUnload <- function(libpath) {
  library.dynam.unload("TreeTools", libpath)
}

## Reminders when releasing for CRAN
release_questions <- function() {
  c(
    "Is the code free of #TODOs?",
    "Have you checked the Vignettes for sanity?"
  )
}


# Additional steps:
#
# Propagate changes in README.md to R/TreeTools-package.R


# Additional tests:
#
# spell_check()
# pkgdown::build_reference_index()
#
# run_examples()
# build_vignettes()
#
# devtools::check_win_devel(quiet = TRUE); rhub::check_for_cran()
# Check valgrind results on Github Actions
# revdepcheck::revdep_check()
#
# codemeta::write_codemeta()
#
# tools::resaveRdaFiles("R", compress="auto") - is default bzip2 the optimal?
# tools::checkRdaFiles("R") - set optimal compression in `data-raw`
