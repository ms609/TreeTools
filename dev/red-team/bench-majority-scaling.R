# Benchmark: majority consensus scaling in k (number of trees) at fixed n.
# The new core scales ~linearly in k; both count modes (hashed default, exact
# opt-in) are timed.
#
# Historical result vs the previous O(k^2 n) core (random trees, n = 100):
#   k:        100    200    400    800   1600
#   new (s):  0.014  0.027  0.050  0.100  0.209   (exponent a = 0.97, linear)
#   old (s):  0.027  0.095  0.353  1.36   5.36    (exponent a = 1.91, quadratic)
#   speed-up: 1.9x   3.6x   7.0x   13.6x  25.6x
#
# Run:  Rscript dev/red-team/bench-majority-scaling.R
suppressMessages(devtools::load_all(".", quiet = TRUE))

# Random trees: their majority consensus is ~unresolved, so the legacy core
# never hits its perfectly-resolved early exit and pays the full O(k^2 n).
n_tip <- 100L
ks <- c(100, 200, 400, 800, 1600)
reps <- 2L
set.seed(1)
make_forest <- function(k) {
  forest <- lapply(seq_len(k), function(i) ape::rtree(n_tip, br = NULL))
  forest <- Preorder(RenumberTips(forest, forest[[1]]))
  structure(forest, class = "multiPhylo")
}

bench1 <- function(fn, forest, p) {
  best <- Inf
  for (i in seq_len(reps)) {
    t0 <- Sys.time()
    fn(forest, p)
    best <- min(best, as.numeric(Sys.time() - t0, units = "secs"))
  }
  best
}

hashed <- function(f, p) TreeTools:::consensus_tree(f, p)
exact  <- function(f, p) TreeTools:::consensus_tree(f, p, exact = TRUE)

cat(sprintf("Majority consensus, fixed n = %d tips, p = 0.5\n", n_tip))
cat(sprintf("%6s %12s %12s\n", "k", "hashed (s)", "exact (s)"))
res <- data.frame()
for (k in ks) {
  forest <- make_forest(k)
  th <- bench1(hashed, forest, 0.5)
  te <- bench1(exact, forest, 0.5)
  cat(sprintf("%6d %12.4f %12.4f\n", k, th, te))
  res <- rbind(res, data.frame(k = k, hashed = th, exact = te))
}

# Scaling exponent: slope of log(time) vs log(k).
fit_exp <- function(k, t) {
  ok <- t > 0
  unname(coef(lm(log(t[ok]) ~ log(k[ok])))[2])
}
cat(sprintf("\nEmpirical scaling exponent (time ~ k^a; target ~1.0, linear):\n"))
cat(sprintf("  hashed a = %.2f\n", fit_exp(res$k, res$hashed)))
cat(sprintf("  exact  a = %.2f\n", fit_exp(res$k, res$exact)))
