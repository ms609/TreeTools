source("benchmark/_init.R") # sets seed

if (!file.exists("benchmark/tr80.rds")) {
  set.seed(1337)
  tr80 <- rtree(80)
  tr2000 <- rtree(2000)
  saveRDS(tr80, "benchmark/tr80.rds")
  saveRDS(tr2000, "benchmark/tr2000.rds")
}

tr80 <- readRDS("benchmark/tr80.rds")
tr2000 <- readRDS("benchmark/tr2000.rds")

rooted80 <- RootTree(tr80, 1)
rooted2k <- RootTree(tr2000, 1)
Benchmark(MatchNodes(tr80, rooted80), times = 500)
Benchmark(MatchNodes(tr2000, rooted2k), times = 50)
Benchmark(MatchEdges(tr80, rooted80), times = 500)

# Preorder() must relocate node labels when nodes are not numbered in preorder
.Scramble <- function(tr) {
  nTip <- NTip(tr)
  nodeOrder <- c(seq_len(nTip + 1), nTip + 1 + sample(tr[["Nnode"]] - 1))
  tr[["edge"]][] <- nodeOrder[tr[["edge"]]]
  attr(tr, "order") <- NULL
  tr[["node.label"]] <- paste0("n", seq_len(tr[["Nnode"]]))
  tr
}
labelled80 <- .Scramble(tr80)
labelled2k <- .Scramble(tr2000)
Benchmark(Preorder(labelled80), times = 500)
Benchmark(Preorder(labelled2k), times = 50)
