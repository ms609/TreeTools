source("benchmark/_init.R")

bigTrees <- as.phylo(1:10, 64*32)
someTrees <- as.phylo(1:240, 64*3)

Benchmark("as.Splits192", ub(as.Splits(someTrees))) # 16.5
Benchmark("as.Splits64.32", ub(as.Splits(bigTrees))) # 23.4

bigSplits <- as.Splits(bigTrees)
someSplits <- as.Splits(someTrees)

Benchmark("as.phylo.some", ub(lapply(someSplits, as.phylo))) # 16.7
Benchmark("as.phylo.big", ub(lapply(bigSplits, as.phylo))) # 34.0
