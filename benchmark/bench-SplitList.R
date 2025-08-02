source("benchmark/_init.R")

bigTrees <- as.phylo(1:10, 64*32)
someTrees <- as.phylo(1:240, 64*3)

Benchmark("as.Splits192", ub(as.Splits(someTrees), times = 100)) # 16.5
Benchmark("as.Splits64.32", ub(as.Splits(bigTrees), times = 100)) # 23.4

bigSplits <- as.Splits(bigTrees)
someSplits <- as.Splits(someTrees)

# Test: -3 or keep power of 2 for SL_MAX_SPLITS?
# See: #define SL_MAX_SPLITS (SL_MAX_TIPS - 3) in SplitList.h
Benchmark("as.phylo.some", ub(lapply(someSplits, as.phylo), times = 180)) # 16.7
Benchmark("as.phylo.big", ub(lapply(bigSplits, as.phylo), times = 180)) # 34.0
