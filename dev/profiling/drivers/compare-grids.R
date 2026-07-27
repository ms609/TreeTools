# Compare two benchmark-grid CSVs (baseline vs optimised) and print per-cell
# speedups. Run after consensus-grid.R has overwritten results-grid.csv.
# Run:  Rscript dev/profiling/drivers/compare-grids.R [old.csv] [new.csv]
a <- commandArgs(TRUE)
oldf <- if (length(a) >= 1) a[1] else "dev/profiling/results-grid-baseline.csv"
newf <- if (length(a) >= 2) a[2] else "dev/profiling/results-grid.csv"
o <- read.csv(oldf); n <- read.csv(newf)
key <- function(d) paste(d$regime, d$scenario, d$n, d$k, d$method, sep = "|")
o$key <- key(o); n$key <- key(n)
m <- merge(o[, c("key", "min_s", "nsplits")],
           n[, c("key", "min_s", "nsplits")], by = "key",
           suffixes = c("_old", "_new"))
m$speedup <- round(m$min_s_old / m$min_s_new, 2)
m$changed_nsplits <- m$nsplits_old != m$nsplits_new   # MUST be all FALSE (correctness)
parts <- do.call(rbind, strsplit(m$key, "|", fixed = TRUE))
m$regime <- parts[, 1]; m$scenario <- parts[, 2]; m$n <- as.integer(parts[, 3])
m$k <- as.integer(parts[, 4]); m$method <- parts[, 5]
m <- m[order(m$method, m$regime, m$n, m$k, m$scenario), ]

cat("== per-cell speedup (old/new min_s); >1 = faster ==\n")
for (i in seq_len(nrow(m))) {
  cat(sprintf("%-7s %-11s %-9s n=%-5d k=%-5d  %7.4f -> %7.4f  x%-5.2f %s\n",
              m$method[i], m$regime[i], m$scenario[i], m$n[i], m$k[i],
              m$min_s_old[i], m$min_s_new[i], m$speedup[i],
              if (m$changed_nsplits[i]) "  *** NSPLITS CHANGED ***" else ""))
}
if (any(m$changed_nsplits)) {
  cat("\n!!! CORRECTNESS REGRESSION: nsplits changed in some cell !!!\n")
} else {
  cat("\nnsplits identical in every cell (no correctness change).\n")
}
cat(sprintf("hashed median speedup: x%.2f; exact median speedup: x%.2f\n",
            median(m$speedup[m$method == "hashed"]),
            median(m$speedup[m$method == "exact"])))
