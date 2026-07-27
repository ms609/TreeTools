# C-004 regression: full Consensus() output must be byte-identical before/after
# replacing the metadata-strip lapply with bare c(trees).
# Run once on the OLD code (saves ref), once on the NEW code (compares).
suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(42)
ref_path <- "dev/profiling/c004-ref.rds"

base  <- lapply(1:9, function(i) ape::rtree(12, br = NULL))
mp    <- structure(RenumberTips(structure(base, class = "multiPhylo"), base[[1]]),
                   class = "multiPhylo")                       # aligned order
mpU   <- structure(base, class = "multiPhylo")                 # UNaligned order
addEL <- function(f) structure(lapply(f, function(t){t$edge.length<-seq_len(nrow(t$edge)); t}), class="multiPhylo")
addNL <- function(f) structure(lapply(f, function(t){t$node.label<-paste0("n",seq_len(t$Nnode)); t}), class="multiPhylo")
diff_sz <- structure(c(base[1:5], lapply(base[6:9], function(t) DropTip(t, 1:2))), class="multiPhylo")

variants <- list(
  plain     = mp,
  unaligned = mpU,                              # exercises check.labels=TRUE relabel
  edge_len  = addEL(mp),
  node_lab  = addNL(mp),
  both      = addNL(addEL(mp)),
  labelled  = ape::.compressTipLabel(mp),
  lab_el    = ape::.compressTipLabel(addEL(mp)),
  list_in   = unclass(mp),                      # plain list, not multiPhylo
  diff_size = diff_sz                           # exercises the KeepTip repeat loop
)
grid <- expand.grid(v = names(variants), p = c(0.5, 2/3, 1),
                    check = c(TRUE, FALSE), hash = c(TRUE, FALSE),
                    stringsAsFactors = FALSE)
canon <- function(tr) {
  if (!inherits(tr, "phylo")) return(tr)
  list(edge = Preorder(RootTree(tr, 1))$edge, tip = tr$tip.label, Nnode = tr$Nnode)
}
out <- Map(function(v, p, check, hash)
  tryCatch(canon(suppressWarnings(Consensus(variants[[v]], p, check, hash))),
           error = function(e) paste("ERR:", conditionMessage(e))),
  grid$v, grid$p, grid$check, grid$hash)

if (!file.exists(ref_path)) {
  saveRDS(out, ref_path)
  cat(sprintf("Saved reference: %d cells -> %s\n", length(out), ref_path))
} else {
  ref <- readRDS(ref_path)
  same <- mapply(identical, ref, out)
  if (all(same)) {
    cat(sprintf("PASS: all %d Consensus() outputs identical to reference.\n", length(out)))
  } else {
    cat(sprintf("FAIL: %d/%d cells differ:\n", sum(!same), length(out)))
    bad <- which(!same)
    for (i in bad) cat(sprintf("  v=%-9s p=%.3g check=%s hash=%s\n",
                               grid$v[i], grid$p[i], grid$check[i], grid$hash[i]))
  }
}
