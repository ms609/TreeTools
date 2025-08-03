source("benchmark/_init.R") # sets seed

tree <- rtree(100)
Benchmark("DropTip", ub(DropTip(tree, 5)))

unif <- tree
unif$edge.length <- NULL
Benchmark("DropTipUnlen", ub(DropTip(unif, 5)))
