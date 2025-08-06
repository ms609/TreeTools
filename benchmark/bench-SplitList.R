source("benchmark/_init.R")

bigTrees <- as.phylo(1:10, 64*32)
someTrees <- as.phylo(1:240, 64*3)

Benchmark(as.Splits(someTrees))
Benchmark(as.Splits(bigTrees))

bigSplits <- as.Splits(bigTrees)
someSplits <- as.Splits(someTrees)

Benchmark(lapply(someSplits, as.phylo))
Benchmark(lapply(bigSplits, as.phylo))
