#' Reweight phylogenetic characters
#' Reweight characters in a phylogenetic dataset
#' 
#' @param dataset A phylogenetic data matrix of \pkg{phangorn} class `phyDat`,
#' or as a matrix in the format produced by [`PhyDatToMatrix()`].
#' @param weights Unnamed integer vector specifying desired weight of each
#' character in turn; or named integer vector specifying weights of each
#' character; unnamed entries will be assigned weight 1.
#' 
#' @return `Reweight()` returns `dataset` after reweighting characters.
#' For a matrix, this will be attained by repeating each column the specified
#' number of times. For a `phyDat` object, the `weights` attribute will be
#' modified.
#' 
#' @examples
#' mat <- rbind(a = c(0, 2, 0), b = c(0, 2, 0), c = c(1, 3, 0), d = c(1, 3, 0))
#' dat <- MatrixToPhyDat(mat)
#' 
#'  # Set character 1 to weight 1, character 2 to weight 2; omit character 3
#' Reweight(mat, c(1, 2, 0))
#' # Equivalently:
#' Reweight(dat, c("3" = 0, "2" = 2))
#' @template MRS
#' @export
Reweight <- function(dataset, weights) {
  UseMethod("Reweight")
}

#' @export
Reweight.matrix <- function(dataset, weights) {
  nChar <- dim(dataset)[[2]]
  
  keys <- as.integer(names(weights))
  if (!is.null(keys) && length(keys) > 0) {
    w <- rep(1, nChar)
    w[as.integer(keys)] <- weights
    weights <- w
  }
  
  if (!is.numeric(weights)) {
    stop("`weights` must be a numeric vector")
  }
  
  if (length(weights) != nChar) {
    stop("Length of `weights` (", length(weights),
         ") must match number of characters (", nChar, ")")
  }
  charName <- colnames(dataset)
  dataset <- dataset[, rep(seq_len(nChar), times = weights), drop = FALSE]
  `colnames<-`(dataset,
               paste0(charName, "_", 
                      unlist(sapply(weights, seq_len, USE.NAMES = FALSE),
                             recursive = FALSE, use.names = FALSE)))
  
}

#' @export
Reweight.phyDat <- function(dataset, weights) {
  
}
