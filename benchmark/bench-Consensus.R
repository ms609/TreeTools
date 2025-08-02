library("TreeDist")

forest1 <- as.phylo(0:200, 80)
forest2 <- as.phylo(0:20, 260)

Benchmark("consensus1", ub(Consensus(forest1)))
Benchmark("consensus2", ub(Consensus(forest2)))

# Benchmark("RF1", ub(RobinsonFoulds(forest1), times = 42))
# Benchmark("RF2", ub(RobinsonFoulds(forest2)))
