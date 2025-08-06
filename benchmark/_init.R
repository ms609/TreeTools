library("TreeTools")
ub <- bench::mark

set.seed(1337)

Benchmark <- function(id, result) {
  if (interactive()) {
    print(result)
  } else {
    saveRDS(result, paste0(id, ".bench.Rds"))
  }
}
