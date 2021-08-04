devtools::load_all()
library("microbenchmark")
set.seed(0)

trees <- Postorder(as.phylo(1:200, 1001))

message(Sys.time(), ": Starting.")



message(Sys.time(), ": End.")
