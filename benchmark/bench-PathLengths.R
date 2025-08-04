source("benchmark/_init.R") # sets seed

tr80 <- rtree(80)
Benchmark("PathLength.r80", ub(PathLengths(tr80))) # 537 µs

tr80Unif <- tr80
tr80Unif[["edge.length"]] <- NULL
Benchmark("PathLength.r80.unif", ub(PathLengths(tr80Unif))) # 522 µs

tr2000 <- rtree(2000)
Benchmark("PathLength.r2000", ub(PathLengths(tr2000), times = 12)) # 461 ms
