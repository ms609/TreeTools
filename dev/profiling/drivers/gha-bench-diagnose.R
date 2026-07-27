# Why do the GHA bench-Consensus.R cells show NSD after C-001/C-002?
# Diagnose: which path each benchmark exercises + wrapper cost breakdown.
suppressMessages(devtools::load_all(".", quiet = TRUE))

forest201.80 <- as.phylo(0:200, 80)
forest21.260 <- as.phylo(0:20, 260)
forest1k.888 <- as.phylo(0:1000, 888)
forestMaj <- c(as.phylo(rep(c(0, 1e9), c(15, 50)), 100),
               structure(lapply(rep(100, 50), BalancedTree), class = "multiPhylo"))

describe <- function(nm, f) {
  cls <- class(f)
  tl  <- attr(f, "TipLabel")
  has_el <- any(vapply(f, function(t) !is.null(t[["edge.length"]]), logical(1)))
  has_nl <- any(vapply(f, function(t) !is.null(t[["node.label"]]), logical(1)))
  # do all trees share tip.label in identical order?
  if (is.null(tl)) {
    labs <- lapply(f, function(t) t[["tip.label"]])
    same_order <- all(vapply(labs[-1], identical, logical(1), labs[[1]]))
  } else same_order <- NA  # labelled: order shared by construction
  cat(sprintf("%-14s k=%-5d n=%-4d class=%-10s TipLabel=%-5s edge.len=%-5s node.lab=%-5s sameTipOrder=%s\n",
              nm, length(f), NTip(f[[1]]), paste(cls, collapse=","),
              !is.null(tl), has_el, has_nl, same_order))
}
cat("== forest structure ==\n")
describe("forest201.80", forest201.80)
describe("forest21.260", forest21.260)
describe("forest1k.888", forest1k.888)
describe("forestMaj",   forestMaj)

# ---- Wrapper cost breakdown (replicate Consensus() internals) -----------------
bk <- function(label, expr, times = 25) {
  expr <- substitute(expr)
  e <- parent.frame()
  t <- replicate(times, { g <- gc(FALSE); s <- Sys.time(); eval(expr, e); as.numeric(Sys.time() - s, units = "secs") })
  cat(sprintf("  %-26s %8.2f ms (min %6.2f)\n", label, 1000*median(t), 1000*min(t)))
  invisible(median(t))
}

breakdown <- function(nm, trees, p = 1, check.labels = TRUE) {
  cat(sprintf("\n== %s  (p=%g, check.labels=%s) ==\n", nm, p, check.labels))
  hash <- TRUE
  bk("full Consensus()", Consensus(trees, p, check.labels, hash))
  bk("metadata-strip lapply", lapply(c(trees), function(tr) { tr[["edge.length"]] <- NULL; tr[["node.label"]] <- NULL; tr }))
  stripped <- lapply(c(trees), function(tr) { tr[["edge.length"]] <- NULL; tr[["node.label"]] <- NULL; tr })
  bk("NTip(trees)", NTip(stripped))
  if (check.labels) bk("RenumberTips", RenumberTips(stripped, stripped[[1]]))
  rt <- if (check.labels) RenumberTips(stripped, stripped[[1]]) else stripped
  bk("Preorder(trees[[1]])", Preorder(rt[[1]]))
  bk("consensus_tree core", TreeTools:::consensus_tree(rt, p, exact = !hash))
}

breakdown("forest1k.888", forest1k.888, p = 1,   check.labels = FALSE)
breakdown("forest201.80", forest201.80, p = 1,   check.labels = FALSE)
breakdown("forest21.260", forest21.260, p = 1,   check.labels = TRUE)
breakdown("forestMaj",    forestMaj,    p = 0.5, check.labels = FALSE)
