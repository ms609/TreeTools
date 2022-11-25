#' @export
as.matrix.phylo <- function(x, ...) {
  sp <- as.Splits(x)
  as.matrix(sp) * 1
}

#' @export
as.logical.phylo <- function(x, ...) {
  as.logical(as.Splits(x))
}

#' @export
as.matrix.Splits <- as.logical.Splits

#' @export
as.matrix.phyDat <- function(x, ambigNA = FALSE, inappNA = ambigNA, ...) {
  PhyDatToMatrix(x, ambigNA, inappNA)
}
