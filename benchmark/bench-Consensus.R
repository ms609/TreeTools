library("TreeDist")

forest1 <- as.phylo(0:200, 80)
forest2 <- as.phylo(0:20, 260)

Benchmark("consensus", Consensus(forest1), Consensus(forest2))
Benchmark("RF", RobinsonFoulds(forest1), RobinsonFoulds(forest2))
