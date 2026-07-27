# Baseline benchmark grid for consensus counting methods.
#
# Times the C++ CORE (TreeTools:::consensus_tree) on forests that are already
# RenumberTips+Preorder'd, so we compare ALGORITHMS apples-to-apples (the shared
# R-wrapper cost is measured separately at the end). Covaries the three regimes
# the user named (high-k/low-n, low-k/high-n, high-k/high-n) with tree SHAPE
# (random vs tall caterpillar) and CONCORDANCE (random / 70%-concordant /
# identical) — the covariates that drive exact's O(k.n.h) and (later) Jansson.
#
# Method dispatch is a single place (`run_method`) so Jansson slots in later.
# Writes a tidy CSV; guards pathological exact cells with a warmup probe.
#
# Run:  Rscript dev/profiling/drivers/consensus-grid.R
suppressWarnings(suppressMessages({
  libpath <- readLines("dev/profiling/.libpath.txt", warn = FALSE)[1]
  library(TreeTools, lib.loc = libpath)
}))
set.seed(5813)

OUT  <- "dev/profiling/results-grid.csv"
REPS <- 4L            # timed reps after 1 warmup
PROBE_BUDGET_S <- 8   # if the warmup rep exceeds this, record it alone & skip reps

timeit <- function(expr) {
  expr <- substitute(expr)
  env  <- parent.frame()
  eval(expr, env)                                   # warmup (also a probe)
  t0 <- Sys.time(); eval(expr, env)
  warm <- as.numeric(Sys.time() - t0, units = "secs")
  if (warm > PROBE_BUDGET_S) return(list(min = warm, med = warm, reps = 1L))
  ts <- numeric(REPS)
  for (i in seq_len(REPS)) {
    t0 <- Sys.time(); eval(expr, env)
    ts[i] <- as.numeric(Sys.time() - t0, units = "secs")
  }
  list(min = min(ts), med = stats::median(ts), reps = REPS)
}

# ---- forest generators (each returns a RenumberTips+Preorder'd multiPhylo) ----
prep <- function(forest) {
  forest <- structure(forest, class = "multiPhylo")
  Preorder(RenumberTips(forest, forest[[1]]))
}
gen <- list(
  rand = function(n, k) prep(lapply(seq_len(k), function(i) ape::rtree(n, br = NULL))),
  tall = function(n, k) { labs <- paste0("t", seq_len(n))
    prep(lapply(seq_len(k), function(i) PectinateTree(sample(labs)))) },
  concord70 = function(n, k) { base <- ape::rtree(n, br = NULL)
    nb <- ceiling(0.7 * k)
    prep(c(rep(list(base), nb),
           lapply(seq_len(k - nb), function(i) ape::rtree(n, br = NULL,
                    tip.label = base$tip.label)))) },
  identical = function(n, k) { base <- ape::rtree(n, br = NULL)
    prep(rep(list(base), k)) }
)

# ---- method dispatch (Jansson added here later) -------------------------------
run_method <- function(forest, p, method) {
  switch(method,
    hashed = TreeTools:::consensus_tree(forest, p, exact = FALSE),
    exact  = TreeTools:::consensus_tree(forest, p, exact = TRUE),
    # jansson = TreeTools:::consensus_tree(forest, p, method = "jansson"),
    stop("unknown method"))
}
METHODS <- c("hashed", "exact")

# ---- grid: (regime, n, k) -----------------------------------------------------
grid <- rbind(
  data.frame(regime = "baseline",      n = 100,   k = 100),
  data.frame(regime = "highK_lowN",    n = 30,    k = 1000),
  data.frame(regime = "highK_lowN",    n = 30,    k = 5000),
  data.frame(regime = "highK_lowN",    n = 50,    k = 2000),
  data.frame(regime = "lowK_highN",    n = 2000,  k = 10),
  data.frame(regime = "lowK_highN",    n = 5000,  k = 10),
  data.frame(regime = "lowK_highN",    n = 10000, k = 5),
  data.frame(regime = "highK_highN",   n = 1000,  k = 300),
  data.frame(regime = "highK_highN",   n = 1000,  k = 1000),
  data.frame(regime = "highK_highN",   n = 2000,  k = 500)
)
SCENARIOS <- names(gen)
P <- 0.5

rows <- list()
first <- TRUE
for (gi in seq_len(nrow(grid))) {
  n <- grid$n[gi]; k <- grid$k[gi]; regime <- grid$regime[gi]
  for (sc in SCENARIOS) {
    tg0 <- Sys.time()
    forest <- gen[[sc]](n, k)
    gen_s <- as.numeric(Sys.time() - tg0, units = "secs")
    for (m in METHODS) {
      tr <- timeit(run_method(forest, P, m))
      res <- run_method(forest, P, m)
      nsplits <- if (is.null(res)) 0L else nrow(res)
      row <- data.frame(regime = regime, scenario = sc, n = n, k = k,
                        method = m, p = P, min_s = round(tr$min, 5),
                        med_s = round(tr$med, 5), reps = tr$reps,
                        nsplits = nsplits, gen_s = round(gen_s, 3))
      rows[[length(rows) + 1]] <- row
      cat(sprintf("%-12s %-9s n=%-5d k=%-5d %-7s min=%8.4f med=%8.4f reps=%d nsplits=%d\n",
                  regime, sc, n, k, m, tr$min, tr$med, tr$reps, nsplits))
      utils::write.table(row, OUT, sep = ",", row.names = FALSE,
                         col.names = first, append = !first)
      first <- FALSE
    }
  }
}

# ---- R-wrapper overhead: Consensus() (full pipeline) vs core-only -------------
cat("\n== wrapper overhead (end-to-end Consensus vs core consensus_tree) ==\n")
for (cfg in list(c(30, 5000), c(100, 100), c(2000, 10), c(1000, 1000))) {
  n <- cfg[1]; k <- cfg[2]
  forest <- gen$rand(n, k)
  # raw (unprepared) copy to feed Consensus(), so its RenumberTips/Preorder run
  raw <- structure(lapply(seq_len(k), function(i) forest[[i]]), class = "multiPhylo")
  te <- timeit(Consensus(raw, p = 0.5, check.labels = TRUE))
  tc <- timeit(TreeTools:::consensus_tree(forest, 0.5, exact = FALSE))
  cat(sprintf("n=%-5d k=%-5d  Consensus()=%8.4f  core=%8.4f  wrapper=%8.4f (%.0f%%)\n",
              n, k, te$min, tc$min, te$min - tc$min, 100 * (te$min - tc$min) / te$min))
}
cat("\nWrote", OUT, "\n")
