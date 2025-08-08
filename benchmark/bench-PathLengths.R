source("benchmark/_init.R") # sets seed

tr80 <- rtree(80)

Benchmark("PathLength.r80-internal",
          ub(TreeTools:::path_lengths(tr80$edge, tr80$edge.length, FALSE)))
Benchmark("PathLength.r80", ub(PathLengths(tr80, full = TRUE)))

tr80Unif <- tr80
tr80Unif[["edge.length"]] <- NULL
Benchmark("PathLength.r80.unif", ub(PathLengths(tr80Unif, full = TRUE)))

tr2000 <- rtree(2000)
Benchmark("PathLength.r2000", ub(PathLengths(tr2000, full = TRUE), times = 42))
