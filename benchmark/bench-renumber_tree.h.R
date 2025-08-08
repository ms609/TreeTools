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

bal40 <- BalancedTree(40)$edge
dbal40 <- Disorder(bal40)
nbal40 <- NodeNumber(dbal40)
Benchmark(TreeTools:::postorder_order(bal40))
Benchmark(TreeTools:::postorder_order(dbal40))
Benchmark(TreeTools:::postorder_order(nbal40))

bal40k <- BalancedTree(40000)$edge
dbal40k <- Disorder(bal40k)
nbal40k <- NodeNumber(dbal40k)
Benchmark(TreeTools:::postorder_order(bal40k))
Benchmark(TreeTools:::postorder_order(dbal40k))
Benchmark(TreeTools:::postorder_order(nbal40k))

pec40 <- PectinateTree(40)$edge
dpec40 <- Disorder(pec40)
npec40 <- NodeNumber(dpec40)
Benchmark(TreeTools:::postorder_order(pec40))
Benchmark(TreeTools:::postorder_order(dpec40))
Benchmark(TreeTools:::postorder_order(npec40))
pec40k <- PectinateTree(40000)$edge
dpec40k <- Disorder(pec40k)
npec40k <- NodeNumber(dpec40k)
Benchmark(TreeTools:::postorder_order(pec40k))
Benchmark(TreeTools:::postorder_order(dpec40k))
Benchmark(TreeTools:::postorder_order(npec40k))

r80 <- RandomTree(80, root = TRUE)
rnd80 <- r80$edge
drnd80 <- Disorder(rnd80)
nrnd80 <- NodeNumber(drnd80)
Benchmark(TreeTools:::postorder_order(rnd80))
Benchmark(TreeTools:::postorder_order(drnd80))
Benchmark(TreeTools:::postorder_order(nrnd80))

