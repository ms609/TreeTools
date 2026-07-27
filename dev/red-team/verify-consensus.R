# Consistency cross-check for the O(kn) consensus / split-frequency machinery on
# adversarial inputs: the hashed (default) and exact split counts must agree, in
# both Consensus() and SplitFrequency().
#
# Historical note: during development this script also compared the new core
# against a transient legacy O(k^2 n) oracle (`consensus_tree_legacy`); that gate
# passed with 0 failures across all fixtures below before the oracle was removed.
# The independent oracle for ongoing regression is `ape::consensus()`, exercised
# by tests/testthat/test-consensus.R.
#
# Run:  Rscript dev/red-team/verify-consensus.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(1)

# Canonical set of splits from a packed RawMatrix (order-independent).
split_set <- function(m, n_tip) {
  if (is.null(m) || nrow(m) == 0) return(character(0))
  out <- character(nrow(m))
  for (r in seq_len(nrow(m))) {
    bits <- as.logical(rawToBits(as.raw(m[r, ])))[seq_len(n_tip)]
    if (isTRUE(bits[1])) bits <- !bits          # tip 1 on the FALSE side
    side <- which(bits)
    if (length(side) < 2 || length(side) > n_tip - 2) next  # skip trivial
    out[r] <- paste(side, collapse = ",")
  }
  sort(unique(out[out != ""]))
}

fails <- 0L
check <- function(cond, msg) {
  if (!isTRUE(cond)) { cat("  FAIL:", msg, "\n"); fails <<- fails + 1L }
}

ps <- c(0.5, 0.55, 0.6, 2/3, 0.75, 0.8, 0.9, 0.99, 1)

# ---- Random forests across a range of (n, k) ----------------------------------
message("== random forests ==")
# k = 1 is handled by R's Consensus() wrapper (returns the single tree), never
# reaching consensus_tree(); the legacy oracle is not robust to it, so skip.
for (n_tip in c(4, 5, 6, 8, 13, 20, 50)) {
  for (k in c(2, 3, 5, 7, 8, 15, 32, 60)) {
    message(sprintf("  n=%d k=%d", n_tip, k))
    forest <- lapply(seq_len(k), function(i) ape::rtree(n_tip, br = NULL))
    forest <- RenumberTips(forest, forest[[1]])
    forest <- Preorder(forest)
    class(forest) <- "multiPhylo"
    for (p in ps) {
      hashed <- split_set(TreeTools:::consensus_tree(forest, p), n_tip)
      exact  <- split_set(TreeTools:::consensus_tree(forest, p, hash = FALSE), n_tip)
      check(identical(hashed, exact),
            sprintf("consensus hashed!=exact n=%d k=%d p=%.3g", n_tip, k, p))
    }
  }
}

# ---- Adversarial structured fixtures -----------------------------------------
cat("== structured fixtures ==\n")
adversarial <- list(
  all_identical = rep(list(BalancedTree(8)), 7),
  star_disagree = list(ape::read.tree(text = "((a,b),(c,d));"),
                       ape::read.tree(text = "((a,c),(b,d));")),
  single_split  = list(ape::read.tree(text = "((a,b,c,d),(e,f,g));"),
                       ape::read.tree(text = "((a,b,c,d),(e,f,g));")),
  one_off_bal   = c(rep(list(BalancedTree(10)), 3), list(PectinateTree(10))),
  threshold_tie = c(rep(list(BalancedTree(8)), 2), rep(list(PectinateTree(8)), 2)),
  mixed_resoln  = list(BalancedTree(8), CollapseNode(BalancedTree(8), 11:12),
                       PectinateTree(8))
)
for (nm in names(adversarial)) {
  forest <- adversarial[[nm]]
  forest <- RenumberTips(forest, forest[[1]])
  forest <- Preorder(forest)
  class(forest) <- "multiPhylo"
  n_tip <- NTip(forest[[1]])
  for (p in ps) {
    hashed <- split_set(TreeTools:::consensus_tree(forest, p), n_tip)
    exact  <- split_set(TreeTools:::consensus_tree(forest, p, hash = FALSE), n_tip)
    check(identical(hashed, exact), sprintf("%s p=%.3g: hashed!=exact", nm, p))
  }
}

# ---- Hashed vs exact split frequencies (gate 5) ------------------------------
cat("== hashed vs exact SplitFrequency ==\n")
for (n_tip in c(4, 6, 8, 13, 30)) {
  for (k in c(1, 2, 5, 12, 40)) {
    forest <- lapply(seq_len(k), function(i) ape::rtree(n_tip, br = NULL))
    class(forest) <- "multiPhylo"
    se <- SplitFrequency(forest, hash = FALSE)
    sh <- SplitFrequency(forest, hash = TRUE)
    # Compare as (split -> count) maps, order-independent
    key <- function(s) {
      if (length(s) == 0) return(character(0))
      paste(as.character(s), attr(s, "count"))
    }
    check(setequal(key(sh), key(se)),
          sprintf("SplitFrequency n=%d k=%d: hashed != exact", n_tip, k))
  }
}

cat(sprintf("\n==== %s : %d failures ====\n",
            if (fails == 0) "PASS" else "FAIL", fails))
quit(status = if (fails == 0) 0 else 1)
