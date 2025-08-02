source("benchmark/_init.R")
bigSplits <- as.Splits(as.phylo(100, 64*32))
someSplits <- as.Splits(as.phylo(100, 64*3))

# Test: -3 or keep power of 2 for SL_MAX_SPLITS?
# See: #define SL_MAX_SPLITS (SL_MAX_TIPS - 3) in SplitList.h
Benchmark("as.Splits", as.phylo(someSplits), as.phylo(bigSplits))
