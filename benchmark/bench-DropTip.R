source("benchmark/_init.R") # sets seed

if (!file.exists("benchmark/tr80.rds")) {
  set.seed(1337)
  tr80 <- rtree(80)
  tr2000 <- rtree(2000)
  saveRDS(tr80, "benchmark/tr80.rds")
  saveRDS(tr2000, "benchmark/tr2000.rds")
}

tr80 <- readRDS("benchmark/tr80.rds")
tr2000 <- readRDS("benchmark/tr2000.rds")

Benchmark(DropTip(tr80, 5))
Benchmark(DropTip(tr2000, 5))

unlen80 <- tr80
unlen80$edge.length <- NULL
unlen2k <- tr2000
unlen2k$edge.length <- NULL
Benchmark(DropTip(unlen80, 5))
Benchmark(DropTip(unlen2k, 5))
