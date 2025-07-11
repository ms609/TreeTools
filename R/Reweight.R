#' Reweight phylogenetic characters
#' 
#' `Reweight()` allows the weights of specific characters in phylogenetic
#' datasets to be arbitrarily adjusted.
#' 
#' This functionality should be employed with care.
#' The underlying principle of parsimony is that all evolutionary steps are
#' equivalent.
#' Setting different weights to different characters is at odds with that
#' principle, so analysis of a re-weighted matrix using a parsimony-based
#' framework is arguably no longer parsimony analysis; on the most permissive
#' view, the criteria used to determine a weighting scheme will always
#' be arbitrary.
#' 
#' It can be useful to relax the criterion that all evolutionary steps are
#' equivalent -- for example, implied weighting 
#' \insertCite{Goloboff1997}{TreeTools} typically recovers better trees than
#' equal-weights parsimony \insertCite{Smith2019}{TreeTools}.
#' This said, assigning different weights to different characters tacitly
#' imposes a model of evolution that differs from that implicit in equal-weights
#' parsimony.  Whereas probabilistic models can be evaluated by various methods
#' (e.g. fit, marginal likelihood, posterior predictive power), there are no
#' principled methods of comparing different models under a parsimony framework.
#' 
#' As such, `Reweight()` is likely to be useful for a narrow set of uses.
#' Examples may include: 
#'  - informal robustness testing, to explore whether certain characters are 
#'  more or less influential on the resulting tree;
#'  - Imposing constraints on a dataset, by adding each constraint as a column
#'  in a dataset whose weight exceeds the total amount of data.
#' 
#' 
#' @param dataset A phylogenetic data matrix of \pkg{phangorn} class `phyDat`,
#' or as a matrix in the format produced by [`PhyDatToMatrix()`].
#' @param weights Unnamed integer vector specifying desired weight of each
#' character in turn; or named integer vector specifying weights of each
#' character; unnamed entries will be assigned weight 1.
#' 
#' @return `Reweight()` returns `dataset` after adjusting the weights of
#' the specified characters.
#' For a matrix, this is attained by repeating each column the `weights` times.
#' For a `phyDat` object, the "weight" attribute will be modified.
#' 
#' @references \insertAllCited{}
#' @examples
#' mat <- rbind(a = c(0, 2, 0), b = c(0, 2, 0), c = c(1, 3, 0), d = c(1, 3, 0))
#' dat <- MatrixToPhyDat(mat)
#' 
#'  # Set character 1 to weight 1, character 2 to weight 2; omit character 3
#' Reweight(mat, c(1, 2, 0))
#' # Equivalently:
#' Reweight(dat, c("3" = 0, "2" = 2))
#' @template MRS
#' @family phylogenetic matrix conversion functions
#' @export
Reweight <- function(dataset, weights) {
  UseMethod("Reweight")
}

.MakeWeights <- function(weights, nChar) {
  if (!is.numeric(weights)) {
    stop("`weights` must be a numeric vector")
  }
  
  keys <- as.integer(names(weights))
  if (!is.null(keys) && length(keys) > 0) {
    if (any(keys < 1 | keys > nChar)) {
      stop("Indices in `names(weights)` must be between 1 and ", nChar, 
           ". Found ", paste(keys[keys < 1 | keys > nChar], collapse = ", "))
    }
    w <- rep(1, nChar)
    w[as.integer(keys)] <- weights
    weights <- w
  }
  
  if (length(weights) != nChar) {
    stop("Length of `weights` (", length(weights),
         ") must match number of characters (", nChar, ")")
  }
  
  # Return:
  weights
}

#' @export
Reweight.matrix <- function(dataset, weights) {
  nChar <- dim(dataset)[[2]]
  weights <- .MakeWeights(weights, nChar)
  
  dataset <- dataset[, rep(seq_len(nChar), times = weights), drop = FALSE]
  `colnames<-`(dataset,
               paste0(colnames(dataset), "_", 
                      unlist(sapply(weights, seq_len, USE.NAMES = FALSE),
                             recursive = FALSE, use.names = FALSE)))
}

#' @export
Reweight.phyDat <- function(dataset, weights) {
  nChar <- length(attr(dataset, "index"))
  weights <- .MakeWeights(weights, nChar)
  wAtt <- double(attr(dataset, "nr"))
  index <- attr(dataset, "index")
  for (i in seq_along(weights)) {
    wAtt[index[[i]]] <- wAtt[index[[i]]] + weights[[i]]
  }
  attr(dataset, "weight") <- wAtt
  
  # Return:
  dataset
}
