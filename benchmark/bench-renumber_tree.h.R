source("benchmark/_init.R") # sets seed

Disorder <- function(edge) {
  set.seed(1)
  edge <- edge[sample(dim(edge)[[1]]), ]
}

NodeNumber <- function(edge) {
  set.seed(1)
  nTip <- (dim(edge)[[1]] + 2) / 2
  nNode <- nTip - 1
  nodes <- nTip + seq_len(nNode)
  newNumber <- sample(nodes)
  edge[edge > nTip] <- - edge[edge > nTip]
  for (i in seq_len(nTip - 1)) {
    edge[edge == -nodes[[i]]] <- newNumber[[i]]
  }
  edge
}

pec40 <- BalancedTree(40)$edge
dpec40 <- Disorder(pec40)
npec40 <- NodeNumber(dpec40)
Benchmark("postorder-pec40", ub(TreeTools:::postorder_order(pec40)))
Benchmark("postorder-dpec40", ub(TreeTools:::postorder_order(dpec40)))
Benchmark("postorder-npec40", ub(TreeTools:::postorder_order(npec40)))

pec40k <- BalancedTree(40000)$edge
dpec40k <- Disorder(pec40k)
npec40k <- NodeNumber(dpec40k)
Benchmark("postorder-pec40k", ub(TreeTools:::postorder_order(pec40k)))
Benchmark("postorder-dpec40k", ub(TreeTools:::postorder_order(dpec40k)))
Benchmark("postorder-npec40k", ub(TreeTools:::postorder_order(npec40k)))

pec40 <- PectinateTree(40)$edge
dpec40 <- Disorder(pec40)
npec40 <- NodeNumber(dpec40)
Benchmark("postorder-pec40", ub(TreeTools:::postorder_order(pec40)))
Benchmark("postorder-dpec40", ub(TreeTools:::postorder_order(dpec40)))
Benchmark("postorder-npec40", ub(TreeTools:::postorder_order(npec40)))
pec40k <- PectinateTree(40000)$edge
dpec40k <- Disorder(pec40k)
npec40k <- NodeNumber(dpec40k)
Benchmark("postorder-pec40k", ub(TreeTools:::postorder_order(pec40k)))
Benchmark("postorder-dpec40k", ub(TreeTools:::postorder_order(dpec40k)))
Benchmark("postorder-npec40k", ub(TreeTools:::postorder_order(npec40k)))

r80 <- RandomTree(80, root = TRUE)
rnd80 <- r80$edge
drnd80 <- Disorder(rnd80)
nrnd80 <- NodeNumber(drnd80)
Benchmark("postorder-rnd80", ub(TreeTools:::postorder_order(rnd80)))
Benchmark("postorder-drnd80", ub(TreeTools:::postorder_order(drnd80)))
Benchmark("postorder-nrnd80", ub(TreeTools:::postorder_order(nrnd80)))

