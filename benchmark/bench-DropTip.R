source("benchmark/_init.R") # sets seed

tr80 <- rtree(80)
cat("Tree hash 80:", digest::digest(tr80), "\n")
tr2000 <- rtree(2000)
cat("Tree hash 2000:", digest::digest(tr2000), "\n")
Benchmark("DropTip.80", ub(DropTip(tr80, 5)))
Benchmark("DropTip.2000", ub(DropTip(tr2000, 5), times = 25))

unlen80 <- tr80
unlen80$edge.length <- NULL
unlen2k <- tr2000
unlen2k$edge.length <- NULL
Benchmark("DropTipUnlen", ub(DropTip(unlen80, 5)))
Benchmark("DropTipUnlen2k", ub(DropTip(unlen2k, 5)))
