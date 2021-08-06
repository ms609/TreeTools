devtools::load_all()
set.seed(0)

trees1k <- Cladewise(as.phylo(1:200, 1001))
treesRl <- ape::read.tree("C:/research/r/rogue-ms/data-raw/simulations/431/all.bs")

message(Sys.time(), ": Starting.")

times <- if (interactive()) 24 else 2000
microbenchmark::microbenchmark(Preorder(trees1k), # 59-66
                               Postorder(trees1k), # 127-143
                               Preorder(treesRl), # 33-37
                               Postorder(treesRl), # 34-38
                               times = times)

message(Sys.time(), ": End.")
