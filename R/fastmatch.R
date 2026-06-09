# Internal drop-in replacements for fastmatch's fmatch() and %fin%.
#
# fastmatch is listed in Suggests, not Imports, because it cannot be compiled to
# WebAssembly: its src/dummy.c declares R_registerRoutines() with a deliberately
# wrong signature (to silence a CRAN check) that fails to link against webR's
# real signature, breaking every TreeTools-dependent Shinylive / webR app at
# load time. Other compiled dependencies (e.g. stringi) are unaffected.
#
# Dispatch is resolved exactly once, in .onLoad() (see zzz.R). When fastmatch is
# installed, .onLoad() rebinds `.FastMatch` and `%fin%` to fastmatch's own function
# objects, so there is no per-call wrapper or branch: the internal names *are*
# fastmatch::fmatch and fastmatch's operator, preserving its cached-hash speed-up
# at zero shim overhead. When fastmatch is absent -- or disabled for testing via
# options(TreeTools.fastmatch = FALSE) -- the base equivalents defined below are
# used. They return results identical to fastmatch for every input.

# Fallback for fastmatch::fmatch(). fmatch()'s formals and semantics match
# base::match() exactly (nomatch = NA_integer_, incomparables = NULL), so the
# base function is a faithful drop-in. `base::` is explicit because R/match.R
# turns `match` into an S4 generic within this namespace; fmatch()/%fin% only
# ever see atomic vectors, so we want the plain base behaviour with no S4
# dispatch, exactly as fastmatch provided.
.FastMatch <- base::match

# Fallback for fastmatch's %fin%, which it defines as
# `fmatch(x, table, nomatch = 0L) > 0L`; base %in% is the identical
# `match(x, table, nomatch = 0L) > 0L`. `base::match` for the same reason as
# above.
`%fin%` <- function(x, table) {
  base::match(x, table, nomatch = 0L) > 0L
}
