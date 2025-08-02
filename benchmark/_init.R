library("TreeTools")
ub <- microbenchmark::microbenchmark
Benchmark <- function(id, result) {
  if (interactive()) {
    print(result)
  } else {
    saveRDS(result, paste0(id, ".bench.Rds"))
  }
}
