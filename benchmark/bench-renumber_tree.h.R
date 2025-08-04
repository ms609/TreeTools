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
bal40k <- BalancedTree(40000)$edge
Benchmark("postorder-bal40", # 1.4, 16.1, 65
          ub(TreeTools:::postorder_order(bal40),
             TreeTools:::postorder_order(Disorder(bal40)),
             TreeTools:::postorder_order(NodeNumber(Disorder(bal40)))))
Benchmark("postorder-bal40k", # 1.4, 16.1, 65
          ub(TreeTools:::postorder_order(bal40k), times = 1,
             TreeTools:::postorder_order(Disorder(bal40k)),
             TreeTools:::postorder_order(NodeNumber(Disorder(bal40k)))))


pec40 <- BalancedTree(40)$edge
pec40k <- BalancedTree(40000)$edge
Benchmark("postorder-pec40", # 1.4, 16.1, 65
          ub(TreeTools:::postorder_order(pec40),
             TreeTools:::postorder_order(Disorder(pec40)),
             TreeTools:::postorder_order(NodeNumber(Disorder(pec40)))))
Benchmark("postorder-pec40k", # 1.4, 16.1, 65
          ub(TreeTools:::postorder_order(pec40k), times = 10,
             TreeTools:::postorder_order(Disorder(pec40k)),
             TreeTools:::postorder_order(NodeNumber(Disorder(pec40k)))))

rnd80 <- RandomTree(80)$edge

Benchmark("postorder-rnd80", # 1.4, 16.1, 65
          ub(TreeTools:::postorder_order(rnd80),
             TreeTools:::postorder_order(Disorder(rnd80)),
             TreeTools:::postorder_order(NodeNumber(Disorder(rnd80)))))
