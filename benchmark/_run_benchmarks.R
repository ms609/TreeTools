benchFiles <- list.files("benchmark", "^bench\\-.*\\.R$", full.names = TRUE)
for (benchFile in benchFiles) {
  system2("Rscript", c("-e", paste0("source('", benchFile, "')")))
}
