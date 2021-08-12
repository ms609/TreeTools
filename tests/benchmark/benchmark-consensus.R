devtools::dev_mode()
recompile <- FALSE
pkgload::load_all(compile = recompile, export_all = FALSE,
                  attach_testthat	= FALSE)
message("Generating trees")
trs <- as.phylo(0:10000, 8888)
message("Generating consensus")
Consensus(trs, p = 0.5)
if (interactive()) {
  trs <- as.phylo(0:1000, 888)
  microbenchmark::microbenchmark(Consensus(trs), # 434 → 404 → 43!
                                 consensus(trs), # 2000...
                                 times = c(12, 1))
}
