source("benchmark/_init.R")
# library("TreeDist")

forest1 <- as.phylo(0:200, 80)
forest2 <- as.phylo(0:20, 260)
forest3 <- as.phylo(0:1000, 888)

Benchmark("consensus1", ub(Consensus(forest1)))
Benchmark("consensus2", ub(Consensus(forest2)))
Benchmark("consensus3", ub(Consensus(forest3)))
# Benchmark("RF1", ub(RobinsonFoulds(forest1), times = 42))
# Benchmark("RF2", ub(RobinsonFoulds(forest2), times = 250))
devtools::dev_mode()

  microbenchmark::microbenchmark(Consensus(trs), # 434 → 404 → 43!
                                 consensus(trs), # 2000...
                                 times = c(12, 1))
}
