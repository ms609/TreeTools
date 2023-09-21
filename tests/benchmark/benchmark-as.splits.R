devtools::load_all()
ub <- microbenchmark::microbenchmark
set.seed(0)

trees1k <- Cladewise(as.phylo(1:200, 1001))

bench <- summary(ub(Preorder(trees1k[[1]]), as.Splits(trees1k[[1]])))
meds <- bench$median
meds[2] / meds[1] # 22.2 -> 8.8

