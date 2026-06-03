#!/usr/bin/env Rscript
# Settle the p = 0.5 exact-tie question empirically:
#   - what does TreeTools::Consensus drop/keep at exactly 50%?
#   - does it agree with ape::consensus at the tie boundary?
#   - does it agree with the documented "proportion p or more" semantics?
suppressMessages(pkgload::load_all(quiet = TRUE))  # run from the package root

cat("TreeTools loaded\n\n")

splitsOf <- function(tr) {
  s <- as.Splits(tr)
  sort(vapply(as.character(s), identity, character(1)))
}

report <- function(label, trees, p) {
  n <- length(trees)
  trees <- RenumberTips(trees, trees[[1]])
  tt <- Consensus(trees, p = p)
  ttx <- Consensus(trees, p = p, hash = FALSE)
  ap <- ape::consensus(trees, p = p)
  thresh <- floor(n * p) + 1
  cat(sprintf("== %s : n=%d, p=%.3f, thresh=floor(n*p)+1=%d (count >= %d to keep)\n",
              label, n, p, thresh, thresh))
  cat(sprintf("   TreeTools Nnode (hashed) = %d ; (exact) = %d ; ape Nnode = %d\n",
              tt$Nnode, ttx$Nnode, ap$Nnode))
  agreeHashExact <- isTRUE(all.equal(RootTree(tt, 1), RootTree(ttx, 1)))
  agreeApe <- isTRUE(all.equal(RootTree(tt, 1), RootTree(ap, 1)))
  cat(sprintf("   hashed==exact : %s ; TreeTools==ape : %s\n",
              agreeHashExact, agreeApe))
  if (!agreeApe) {
    cat("   TreeTools splits:\n"); print(splitsOf(tt))
    cat("   ape splits:\n");       print(splitsOf(ap))
  }
  cat("\n")
  invisible(list(tt = tt, ape = ap, agreeApe = agreeApe))
}

# 1. The canonical even tie: 2 + 2 conflicting topologies (the test-Consensus tie fixture)
tie <- c(rep(list(BalancedTree(8)), 2L), rep(list(PectinateTree(8)), 2L))
report("tie 2xBal + 2xPec (8 tip)", tie, 0.5)

# 2. Minimal 2-tree conflict (n = 2, each conflicting split at exactly 50%)
two <- list(ape::read.tree(text = "((a, b), (c, d));"),
            ape::read.tree(text = "((a, c), (b, d));"))
report("2 conflicting quartets", two, 0.5)

# 3. A controlled exact-50% on ONE clade: 4 trees, 2 share clade {t1,t2}, 2 share {t2,t3}
#    (rest identical) so exactly one pair of conflicting splits sits at 2/4.
mk <- function(txt) ape::read.tree(text = txt)
ctrl <- list(
  mk("(((t1,t2),t3),t4,t5,t6);"),
  mk("(((t1,t2),t3),t4,t5,t6);"),
  mk("((t1,(t2,t3)),t4,t5,t6);"),
  mk("((t1,(t2,t3)),t4,t5,t6);")
)
report("controlled 2-2 clade tie (6 tip)", ctrl, 0.5)

# 4. Threshold boundary at p = 2/3 with n = 3 (n*p = 2 exactly): is a 2/3 split kept?
#    doc 'p or more' => count/n >= 2/3 => count >= 2 keep; impl floor(2)+1=3 => drop.
thr <- list(
  mk("(((t1,t2),t3),t4,t5,t6);"),
  mk("(((t1,t2),t3),t4,t5,t6);"),
  mk("((t1,(t2,t3)),t4,t5,t6);")
)
report("p=2/3, n=3, clade at exactly 2/3", thr, 2 / 3)

cat("done\n")
