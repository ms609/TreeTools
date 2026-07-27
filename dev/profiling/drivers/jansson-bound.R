# Measured lower bound on Jansson (Maj_Rule_Plus) cost, without implementing it.
#
# Bound (advisor-validated):
#   Phase 2 of Maj_Rule_Plus ~= k Day-matchings of the candidate tree against
#     each input ~= the STRICT path (consensus_tree p=1): k table builds + k O(n)
#     matching passes.
#   Phase 1 >= Phase 2 (same matching PLUS One_Way_Compatible + Merge_Trees +
#     delete/insert + a candidate-table rebuild each iteration).
#   => Jansson total >= 2 * strict-path.
# Strict's inner loop (single-reference contiguity + CLUSTONL/R, no hash map) is
# TIGHTER than a real Phase 2, so this UNDER-counts Jansson => a "2*strict >=
# hashed" verdict is conservative.
#
# `rand` forests have NO strict early-exit (strict consensus is a star: 0 of n-3
# splits found, never hits the perfectly-resolved return) => k full passes, the
# worst case for strict. We also test `tall`.
#
# Verdict per cell: if 2*strict >= hashed, Jansson cannot win there.
# Run:  Rscript dev/profiling/drivers/jansson-bound.R
suppressWarnings(suppressMessages({
  libpath <- readLines("dev/profiling/.libpath.txt", warn = FALSE)[1]
  library(TreeTools, lib.loc = libpath)
}))
set.seed(5813)
REPS <- 4L
timeit <- function(expr) {
  expr <- substitute(expr); env <- parent.frame()
  eval(expr, env)
  ts <- numeric(REPS)
  for (i in seq_len(REPS)) { t0 <- Sys.time(); eval(expr, env)
    ts[i] <- as.numeric(Sys.time() - t0, units = "secs") }
  min(ts)
}
prep <- function(f) { f <- structure(f, class = "multiPhylo"); Preorder(RenumberTips(f, f[[1]])) }
gen <- list(
  rand = function(n, k) prep(lapply(seq_len(k), function(i) ape::rtree(n, br = NULL))),
  tall = function(n, k) { labs <- paste0("t", seq_len(n))
    prep(lapply(seq_len(k), function(i) PectinateTree(sample(labs)))) },
  # 70% identical base + 30% random: strict consensus is still ~star (the random
  # 30% kill the base's splits) so strict does k full passes => bound stays valid,
  # while distinct-split count is bounded => a more realistic-leaning case.
  concord70 = function(n, k) { base <- ape::rtree(n, br = NULL); nb <- ceiling(0.7 * k)
    prep(c(rep(list(base), nb), lapply(seq_len(k - nb), function(i)
      ape::rtree(n, br = NULL, tip.label = base$tip.label)))) },
  # MODERATE conflict (realistic gene-tree / bootstrap-with-rogues): each
  # replicate = base with ~5% of taxa (rogues) relocated to random edges. The
  # majority IS well resolved (rogues each break only a few replicates' backbone
  # splits) but strict stays near-star (some replicate breaks each split) => no
  # early-exit => bound valid. Distinct-split count ~ n + O(k*rogues): moderate.
  moderate = function(n, k) { base <- ape::rtree(n, br = NULL)
    m <- min(25L, max(3L, n %/% 50L))
    prep(lapply(seq_len(k), function(i) {
      rogues <- sample(base[["tip.label"]], m)
      t <- DropTip(base, rogues)
      for (r in rogues) t <- AddTip(t, label = r)  # default where = random edge
      t
    })) }
)
cells <- list(c(30,5000), c(50,2000), c(5000,10), c(10000,5), c(2000,500), c(1000,1000))

cat(sprintf("%-6s %-5s %-6s %10s %10s %12s %s\n",
            "shape","n","k","hashed","strict","2*strict","verdict"))
worst_ratio <- Inf
for (sc in c("rand","tall","concord70","moderate")) for (cell in cells) {
  n <- cell[1]; k <- cell[2]
  f <- gen[[sc]](n, k)
  th <- timeit(TreeTools:::consensus_tree(f, 0.5, exact = FALSE))   # hashed majority
  ts <- timeit(TreeTools:::consensus_tree(f, 1))                    # strict (Jansson lower-bound unit)
  bound <- 2 * ts
  verdict <- if (bound >= th) "Jansson can't win" else "INCONCLUSIVE — implement"
  worst_ratio <- min(worst_ratio, bound / th)
  cat(sprintf("%-6s %-5d %-6d %10.4f %10.4f %12.4f  %s (2*strict/hashed=%.2f)\n",
              sc, n, k, th, ts, bound, verdict, bound / th))
}
cat(sprintf("\nWorst-case 2*strict/hashed across all cells: %.2f\n", worst_ratio))
cat(if (worst_ratio >= 1) "==> CONCLUSIVE: Jansson cannot beat hashed in any tested regime.\n"
    else "==> Some cell inconclusive; implement Jansson there and measure.\n")
