source("benchmark/_init.R")

forest201.80 <- as.phylo(0:200, 80)
forest21.260 <- as.phylo(0:20, 260)
forest1k.888 <- as.phylo(0:1000, 888)
forestMaj <- c(as.phylo(rep(c(0, 1e9), c(15, 50)), 100),
               structure(lapply(rep(100, 50), BalancedTree), class = "multiPhylo"))

Benchmark(Consensus(forest21.260))
Benchmark(Consensus(forest21.260, 0.5, FALSE))
Benchmark(Consensus(forest201.80, check = FALSE))
Benchmark(Consensus(forest1k.888, check = FALSE))

Benchmark(Consensus(forestMaj, 0.5, FALSE))


# library("TreeDist")
# Benchmark("RF1", ub(RobinsonFoulds(forest1), times = 42))
# Benchmark("RF2", ub(RobinsonFoulds(forest2), times = 250))
