ub <- microbenchmark::microbenchmark
pv <- profvis::profvis

devtools::load_all()

bigSplits <- as.Splits(as.phylo(100, 64*32))
someSplits <- as.Splits(as.phylo(100, 64*3))

# Test: -3 or keep power of 2 for SL_MAX_SPLITS?
# See: #define SL_MAX_SPLITS (SL_MAX_TIPS - 3) in SplitList.h
ub(as.phylo(someSplits), as.phylo(bigSplits))
# No replicable difference.
