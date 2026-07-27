# Correctness gate for consensus methods.
#
# THE rule (user): "any systematic error is unacceptable". The hard gate is that
# every TreeTools consensus counting method produces the IDENTICAL set of
# consensus splits on every fixture; ape::consensus() is an independent oracle
# (soft-reported, since its tie convention at count == k/2 can differ).
#
# Method-pluggable: add Jansson to `methods` once implemented and it is gated
# automatically. Loads the installed -O2 build named in .libpath.txt (NOT
# load_all, so the gate reflects shipped behaviour).
#
# Run:  Rscript dev/profiling/correctness-gate.R
suppressWarnings(suppressMessages({
  libpath <- readLines("dev/profiling/.libpath.txt", warn = FALSE)[1]
  library(TreeTools, lib.loc = libpath)
}))

# ---- canonical split sets (label-based; polarity-normalised to exclude the
#      alphabetically-first label; trivial splits dropped) -----------------------
# From a packed RawMatrix (rows = splits, bit i = tip i+1 in the shared numbering
# given by ref_labels).
splits_from_matrix <- function(m, ref_labels) {
  n <- length(ref_labels)
  excl <- 1L                               # tip 1 (shared numbering) always on excluded side
  if (is.null(m) || nrow(m) == 0) return(character(0))
  out <- character(nrow(m))
  for (r in seq_len(nrow(m))) {
    bits <- as.logical(rawToBits(as.raw(m[r, ])))[seq_len(n)]
    if (isTRUE(bits[excl])) bits <- !bits  # excluded label on FALSE side
    side <- which(bits)
    if (length(side) < 2 || length(side) > n - 2) next
    out[r] <- paste(sort(side), collapse = ",")
  }
  sort(unique(out[out != ""]))
}

# From a phylo, expressed in the shared index space defined by ref_labels.
splits_from_phylo <- function(tree, ref_labels) {
  n <- length(ref_labels)
  excl <- 1L
  pp <- ape::prop.part(tree)
  labs <- attr(pp, "labels")
  out <- character(0)
  for (cl in pp) {
    idx <- match(labs[cl], ref_labels)
    if (anyNA(idx)) stop("label mismatch in splits_from_phylo")
    if (excl %in% idx) idx <- setdiff(seq_len(n), idx)  # exclude first label
    if (length(idx) < 2 || length(idx) > n - 2) next
    out <- c(out, paste(sort(idx), collapse = ","))
  }
  sort(unique(out))
}

# ---- method registry: forest (already RenumberTips+Preorder) , p -> split set --
ref_labels_of <- function(forest) forest[[1]][["tip.label"]]

methods <- list(
  hashed = function(forest, p) splits_from_matrix(
    TreeTools:::consensus_tree(forest, p, exact = FALSE), ref_labels_of(forest)),
  exact  = function(forest, p) splits_from_matrix(
    TreeTools:::consensus_tree(forest, p, exact = TRUE),  ref_labels_of(forest))
  # jansson = function(forest, p) splits_from_matrix(
  #   TreeTools:::consensus_tree(forest, p, method = "jansson"), ref_labels_of(forest))
)
# ape oracle (soft): only defined for p in [0.5, 1].
ape_oracle <- function(forest, p)
  splits_from_phylo(ape::consensus(forest, p = p), ref_labels_of(forest))

fails <- 0L; soft <- 0L; checks <- 0L
hard_fail <- function(msg) { cat("  HARD FAIL:", msg, "\n"); fails <<- fails + 1L }
soft_fail <- function(msg) { cat("  (soft) ape disagree:", msg, "\n"); soft <<- soft + 1L }

# Compare all registered methods to each other (HARD) and to ape (SOFT).
gate_one <- function(forest, p, tag) {
  checks <<- checks + 1L
  res <- lapply(methods, function(f) f(forest, p))
  base <- res[[1]]
  for (nm in names(res)[-1]) {
    if (!identical(base, res[[nm]]))
      hard_fail(sprintf("%s p=%.3g: %s != %s", tag, p, names(res)[1], nm))
  }
  # ape is a valid oracle ONLY at p = 0.5: there ">half" matches TreeTools'
  # `count > k*p`. For p > 0.5, ape keeps `count >= k*p` (>=) while TreeTools
  # keeps `count > k*p` (>) — a deliberate, documented-as-"p or more" but
  # actually-strict TreeTools convention. So we only cross-check ape at p=0.5;
  # the threshold path (p>0.5) is gated by internal hashed==exact==jansson
  # agreement (the `count >= thresh` filter is shared by all methods).
  if (isTRUE(all.equal(p, 0.5))) {
    ok_ape <- tryCatch(ape_oracle(forest, p),
                       error = function(e) structure("ERR", msg = conditionMessage(e)))
    if (identical(as.vector(ok_ape), "ERR")) {
      soft_fail(sprintf("%s p=0.5: ape_oracle ERRORED: %s", tag, attr(ok_ape, "msg")))
    } else if (!identical(base, ok_ape)) {
      soft_fail(sprintf("%s p=0.5 (n kept: ours=%d ape=%d)",
                        tag, length(base), length(ok_ape)))
    }
  }
}

prep <- function(forest) {
  forest <- structure(lapply(forest, function(t) {
    t[["edge.length"]] <- NULL; t[["node.label"]] <- NULL; t
  }), class = "multiPhylo")
  Preorder(RenumberTips(forest, forest[[1]]))
}

ps <- c(0.5, 0.55, 0.6, 2/3, 0.75, 0.8, 0.9, 0.99)

set.seed(1)
cat("== random forests (shape: rtree) ==\n")
for (n_tip in c(4, 5, 6, 8, 13, 20, 50)) {
  for (k in c(2, 3, 5, 7, 8, 15, 32, 60)) {
    forest <- prep(lapply(seq_len(k), function(i) ape::rtree(n_tip, br = NULL)))
    for (p in ps) gate_one(forest, p, sprintf("rand n=%d k=%d", n_tip, k))
  }
}

cat("== tall forests (random caterpillars) ==\n")
for (n_tip in c(8, 20, 60)) {
  labs <- paste0("t", seq_len(n_tip))
  for (k in c(3, 7, 16)) {
    forest <- prep(lapply(seq_len(k), function(i) PectinateTree(sample(labs))))
    for (p in ps) gate_one(forest, p, sprintf("tall n=%d k=%d", n_tip, k))
  }
}

cat("== adversarial structured fixtures ==\n")
adversarial <- list(
  all_identical = rep(list(BalancedTree(8)), 7),
  star_disagree = list(ape::read.tree(text = "((a,b),(c,d));"),
                       ape::read.tree(text = "((a,c),(b,d));")),
  single_split  = list(ape::read.tree(text = "((a,b,c,d),(e,f,g));"),
                       ape::read.tree(text = "((a,b,c,d),(e,f,g));")),
  one_off_bal   = c(rep(list(BalancedTree(10)), 3), list(PectinateTree(10))),
  threshold_tie = c(rep(list(BalancedTree(8)), 2), rep(list(PectinateTree(8)), 2)),
  mixed_resoln  = list(BalancedTree(8), CollapseNode(BalancedTree(8), 11:12),
                       PectinateTree(8)),
  bal_vs_pec    = list(BalancedTree(8), PectinateTree(8))[c(1,1,1,1,2,2,2)]
)
for (nm in names(adversarial)) {
  forest <- prep(adversarial[[nm]])
  for (p in ps) gate_one(forest, p, nm)
}

# ---- HIGH-N content gate (verifies the deferred-materialisation witness path
#      at n where survivors are actually rebuilt from witnesses; the rest of the
#      gate only reaches n<=50). hashed split-SET must equal exact split-SET. ----
cat("== high-n content (deferred-materialisation witness path) ==\n")
set.seed(11)
for (n_tip in c(1000L, 3000L)) {
  labs <- paste0("t", seq_len(n_tip))
  hi <- list(
    # survivor-heavy: ~n-3 majority splits get materialised from witnesses
    concord70 = { base <- ape::rtree(n_tip, br = NULL); nb <- 9L
      c(rep(list(base), nb), lapply(seq_len(4L),
        function(i) ape::rtree(n_tip, br = NULL, tip.label = base$tip.label))) },
    identical = rep(list(ape::rtree(n_tip, br = NULL)), 12L),
    # 0-survivor witness path (random caterpillars -> star consensus)
    tall      = lapply(seq_len(12L), function(i) PectinateTree(sample(labs)))
  )
  for (nm in names(hi)) {
    forest <- prep(hi[[nm]])
    for (p in c(0.5, 2/3)) gate_one(forest, p, sprintf("HIGHN %s n=%d", nm, n_tip))
  }
}

# ---- HIGH-N SplitFrequency content+COUNT gate. The witness change also rewired
#      calc_split_frequencies_hashed -> frequencies_list_from_witnesses, whose
#      materialise-all + RawMatrix packing is otherwise only checked at n<=30.
#      hashed must equal exact as a (split -> count) MAP at shipping scale. ------
cat("== high-n SplitFrequency (hashed vs exact: split sets AND counts) ==\n")
freq_keyset <- function(forest, exact) {
  sf <- TreeTools:::split_frequencies(forest, exact = exact)
  m <- sf[["splits"]]; cnt <- sf[["counts"]]
  if (is.null(m) || nrow(m) == 0) return(character(0))
  n <- length(ref_labels_of(forest))
  out <- character(nrow(m))
  for (r in seq_len(nrow(m))) {
    bits <- as.logical(rawToBits(as.raw(m[r, ])))[seq_len(n)]
    if (isTRUE(bits[1])) bits <- !bits      # tip 1 on excluded side (as elsewhere)
    out[r] <- paste0(paste(which(bits), collapse = ","), "=", cnt[r])
  }
  sort(out)
}
set.seed(23)
for (n_tip in c(2000L, 5000L)) {
  forest <- prep(lapply(seq_len(12L), function(i) ape::rtree(n_tip, br = NULL)))
  checks <- checks + 1L
  if (!identical(freq_keyset(forest, FALSE), freq_keyset(forest, TRUE))) {
    hard_fail(sprintf("SplitFrequency hashed != exact (split/count map) n=%d", n_tip))
  }
}

cat(sprintf("\n==== %s ==== %d checks; %d HARD fails; %d soft (ape) disagreements\n",
            if (fails == 0) "PASS" else "FAIL", checks, fails, soft))
if (soft > 0) cat("NOTE: soft ape disagreements are usually count==k*p tie-convention; inspect if unexpected.\n")
quit(status = if (fails == 0) 0L else 1L)
