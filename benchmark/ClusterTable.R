ub <- microbenchmark::microbenchmark

devtools::load_all()

set.seed(100)
t10 <- rtree(10)
t100 <- rtree(100)
t1k <- rtree(1000)


ub(as.ClusterTable(t10),
   as.ClusterTable(t100),
   as.ClusterTable(t1k)) # edge lengths = 10x slower

# expr    min      lq     mean  median      uq    max neval
# as.ClusterTable(t10)   61.7   66.20   71.298   70.20   73.35  137.3   100
# as.ClusterTable(t100)  164.6  171.25  177.553  174.50  179.85  242.9   100
# as.ClusterTable(t1k) 1117.3 1133.75 1154.081 1150.85 1171.75 1229.8   100

set.seed(444)
trees0 <- lapply(seq_len(10), function(x) rtree(10))
trees1 <- lapply(seq_len(10), function(x) rtree(100))
trees2 <- lapply(seq_len(10), function(x) rtree(1000))
ub(TreeDist::RobinsonFoulds(trees0), times = 1000) # 1.67 with changes
ub(TreeDist::RobinsonFoulds(trees1), TreeDist::RobinsonFoulds(trees2), times = 1000)


#                            expr    min     lq     mean median      uq    max neval
# TreeDist::RobinsonFoulds(trees) 2.7682 2.8622 2.950523 2.9699 3.01405 3.4876   100
# median 2.71 with changes
# Originally 2.8, 12.8
