source("benchmark/_init.R") # sets seed

tr80 <- rtree(80)

Benchmark(TreeTools:::path_lengths(tr80$edge, tr80$edge.length, FALSE),
          min_time = 2)
Benchmark(PathLengths(tr80, full = TRUE), min_time = 2)

tr80Unif <- tr80
tr80Unif[["edge.length"]] <- NULL
Benchmark(PathLengths(tr80Unif, full = TRUE), min_time = 2)

tr2000 <- rtree(2000)
Benchmark(PathLengths(tr2000, full = TRUE))
