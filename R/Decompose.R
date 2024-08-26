#' Decompose additive (ordered) phylogenetic characters
#' 
#' `Decompose()` decomposes additive characters into a series of binary
#' characters, which is mathematically equivalent when analysed under
#' equal weights parsimony.  (This equivalence is not exact
#' under implied weights or under probabilistic tree inference methods.)
#' 
#' @template datasetParam
#' @param indices Integer or logical vector specifying indices of characters
#' that should be decomposed
#' 
#' @return `Decompose()` returns a `phyDat` object in which the specified
#' characters have been decomposed into a number of binary characters.  
#' 
#' @examples
#' 
#' 
#' @template MRS
#' @export
Decompose <- function(dataset, indices) {
  
}
