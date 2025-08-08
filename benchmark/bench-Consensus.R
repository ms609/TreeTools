source("benchmark/_init.R")

forest201.80 <- as.phylo(0:200, 80)
forest21.260 <- as.phylo(0:20, 260)
forest1k.888 <- as.phylo(0:1000, 888)

Benchmark(Consensus(forest21.260))
Benchmark(Consensus(forest201.80))
Benchmark(Consensus(forest1k.888))


# library("TreeDist")
# Benchmark("RF1", ub(RobinsonFoulds(forest1), times = 42))
# Benchmark("RF2", ub(RobinsonFoulds(forest2), times = 250))
