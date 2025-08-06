source("benchmark/_init.R")

forest1 <- as.phylo(0:200, 80)
forest2 <- as.phylo(0:20, 260)
forest3 <- as.phylo(0:1000, 888)

Benchmark(Consensus(forest1))
Benchmark(Consensus(forest2))
Benchmark(Consensus(forest3))


# library("TreeDist")
# Benchmark("RF1", ub(RobinsonFoulds(forest1), times = 42))
# Benchmark("RF2", ub(RobinsonFoulds(forest2), times = 250))
