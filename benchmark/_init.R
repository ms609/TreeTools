library("TreeTools")
Benchmark <- function(id, ...) {
  results <- microbenchmark::microbenchmark(...)
  if (interactive()) {
    print(results)
  } else {
    saveRDS(results, paste0(id, ".bench.Rds"))
  }
}
