devtools::load_all()
set.seed(0)

trees <- Postorder(as.phylo(1:200, 1001))

message(Sys.time(), ": Starting.")

times <- if (interactive()) 10 else 1000
microbenchmark::microbenchmark(Preorder(trees), times = times) # 204-237

message(Sys.time(), ": End.")
