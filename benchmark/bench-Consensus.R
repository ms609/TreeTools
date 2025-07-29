# Hoping for gains when refactoring ClusterTable.h
library("TreeTools")
library("TreeDist")

forest1 <- as.phylo(0:200, 80)
forest2 <- as.phylo(0:20, 260)

ub(Consensus(forest1), Consensus(forest2),
   RobinsonFoulds(forest1), RobinsonFoulds(forest2))

p <- 0.5
benchmark_result <- benchmark_consensus_simple(trees, p, n_iterations = 5)
print(benchmark_result$results_match)  # Sh

# INSTALLED! With Xarr = Rcpp::IntegerMatrix
Unit: milliseconds
                   expr     min       lq      mean  median       uq     max neval
     Consensus(forest1)  5.8774  6.07130  6.767076  6.1733  6.38860 24.4763   100
     Consensus(forest2)  1.8291  1.92865  2.143908  2.0701  2.19085  7.0781   100
RobinsonFoulds(forest1) 44.8181 45.63170 47.383893 46.0655 47.79765 55.4328   100
RobinsonFoulds(forest2)  3.9119  3.99055  4.229728  4.0445  4.10370 10.0436   100

# LOAD_ALL: With Xarr = Rcpp::IntegerMatrix
Unit: milliseconds
                   expr     min       lq      mean   median       uq     max neval
     Consensus(forest1) 14.1216 14.39315 15.403004 14.52975 14.72175 54.8197   100
     Consensus(forest2)  7.2951  7.36435  7.669128  7.41025  7.47750 18.2175   100
RobinsonFoulds(forest1) 44.5817 45.39890 47.318351 45.69750 46.62185 63.7717   100
RobinsonFoulds(forest2)  3.8466  3.97205  4.322742  4.01715  4.10060 11.0708   100
