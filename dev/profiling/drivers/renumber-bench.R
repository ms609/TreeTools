# Benchmark the RenumberTips unlabelled fast path (C++ renumber_tips_to) against
# the per-tree R path it replaces, and the resulting Consensus() wrapper share.
# Run:  Rscript dev/profiling/drivers/renumber-bench.R
suppressWarnings(suppressMessages({
  libpath <- readLines("dev/profiling/.libpath.txt", warn = FALSE)[1]
  library(TreeTools, lib.loc = libpath)
}))
set.seed(5813)
REPS <- 6L
tt <- function(expr) { e <- substitute(expr); env <- parent.frame()
  eval(e, env); ts <- numeric(REPS)
  for (i in seq_len(REPS)) { t0 <- Sys.time(); eval(e, env)
    ts[i] <- as.numeric(Sys.time() - t0, units = "secs") }
  min(ts) }

# Per-tree R reference path (what RenumberTips.multiPhylo used to do for an
# unlabelled forest): apply RenumberTips to each tree individually.
r_path <- function(forest, target) lapply(forest, RenumberTips, target)

cat(sprintf("%-6s %-6s %12s %12s %9s\n", "n", "k", "R per-tree", "C++ batch", "speedup"))
for (cfg in list(c(30,5000), c(50,2000), c(100,1000), c(500,300), c(2000,100))) {
  n <- cfg[1]; k <- cfg[2]
  forest <- lapply(seq_len(k), function(i) ape::rtree(n, br = NULL))  # shuffled orders
  target <- paste0("t", seq_len(n))
  mp <- structure(forest, class = "multiPhylo")
  tR <- tt(r_path(forest, target))
  tC <- tt(RenumberTips(mp, target))
  cat(sprintf("%-6d %-6d %12.4f %12.4f %8.2fx\n", n, k, tR, tC, tR / tC))
}

cat("\n== Consensus() end-to-end vs core (new wrapper share) ==\n")
for (cfg in list(c(30,5000), c(100,100), c(2000,10))) {
  n <- cfg[1]; k <- cfg[2]
  forest <- lapply(seq_len(k), function(i) ape::rtree(n, br = NULL))
  raw <- structure(forest, class = "multiPhylo")
  pre <- Preorder(RenumberTips(raw, raw[[1]]))
  te <- tt(Consensus(raw, p = 0.5))
  tc <- tt(TreeTools:::consensus_tree(pre, 0.5, exact = FALSE))
  cat(sprintf("n=%-5d k=%-5d  Consensus()=%8.4f  core=%8.4f  wrapper=%8.4f (%.0f%%)\n",
              n, k, te, tc, te - tc, 100 * (te - tc) / te))
}
