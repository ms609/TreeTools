source("benchmark/_init.R")

forest201.80 <- as.phylo(0:200, 80)
forest21.888 <- as.phylo(0:20, 888)

Benchmark(TreeDist::RobinsonFoulds(forest201.80))
Benchmark(TreeDist::RobinsonFoulds(forest21.888))
