#' Ratio of external:internal edge length
#' 
#' Reports the ratio of tree length associated with external edges (i.e.
#' edges whose child is a leaf) and internal edges.
#' Where tree length is dominated by internal edges, variation between tips
#' is dominantly controlled by phylogenetic history.
#' 
#' @param x A tree of class \code{\link[ape:read.tree]{phylo}}.
#' @returns `EdgeRatio()` returns a numeric specifying the ratio of external
#' to internal edge length (> 1 means the length of a tree is predominantly
#' in external edges), with atttributes `external`, `internal`, and `total`
#' speficying the total length associated with edges of that nature.
#' @template MRS
#' @export
EdgeRatio <- function(x) UseMethod("EdgeRatio")
  
#' @rdname EdgeRatio
#' @export
EdgeRatio.phylo <- function(x) {
  el <- x[["edge.length"]]
  if (is.null(el)) {
    warning("Edge lengths not specified")
    return(NA_real_)
  }
  ed <- x[["edge"]]
  nTip <- NTip(x)
  external <- ed[, 2] <= nTip
  exLen <- sum(el[external])
  inLen <- sum(el[!external])
  structure(exLen / inLen,
            external = exLen,
            internal = inLen,
            total = sum(exLen, inLen))
}
