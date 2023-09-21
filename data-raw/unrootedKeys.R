suppressPackageStartupMessages(devtools::load_all())

.UnrootedKeys <- function(nTip) {
  if (nTip > 28L) {
    stop("Too many shapes to calculate with ", nTip, " tips.")
  } else if (nTip > 5L) {
    #TODO make efficient - this is horrible!
    shapes <- as.integer(structure(
      vapply(seq_len(as.integer(NRootedShapes(nTip))) - 1L,
             function(shape) UnrootedTreeKey(RootedTreeWithShape(shape, nTip)),
             integer64(1L)),
      class = "integer64"))
    uniqueShapes <- unique(shapes)
  } else {
    uniqueShapes <- 0
  }

  # Return:
  sort(uniqueShapes)
}

data("unrootedKeys", package = "TreeTools")
message("Calculating unrootedKeys for ", length(unrootedKeys) + 1L, " tips.")
unrootedKeys <- c(unrootedKeys, list(.UnrootedKeys(length(unrootedKeys) + 1L)))
message(length(unrootedKeys[[length(unrootedKeys)]]), " keys found.")

usethis::use_data(unrootedKeys, overwrite = TRUE, compress = "xz")
