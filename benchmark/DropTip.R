ub <- microbenchmark::microbenchmark
pv <- profvis::profvis

devtools::load_all()

set.seed(100)
tree <- rtree(100)
unif <- tree
unif$edge.length <- NULL

ub(DropTip(tree, 5), DropTip(unif, 5)) # edge lengths = 10x slower
pv(replicate(10000, DropTip.phylo(tree, 5))) # Bottleneck: path_lengths
