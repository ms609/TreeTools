devtools::load_all()
set.seed(0)

trees1k <- Cladewise(as.phylo(1:200, 1001))
treesRl <- ape::read.tree("C:/research/r/rogue-ms/data-raw/simulations/431/all.bs")

message(Sys.time(), ": Starting.")

mb <- microbenchmark::microbenchmark
times <- if (interactive()) 24 else 200
mb(Preorder(trees1k), # 59-66 --> 38
                               Postorder(trees1k), # 127-143 --> 51
                               Preorder(treesRl), # 33-37 --> 23
                               Postorder(treesRl), # 34-38 --> 30
                               times = times)

ub(Postorder(trees1k[[1]]$edge),
   postorder_order(trees1k[[1]]$edge))

message(Sys.time(), ": End.")
