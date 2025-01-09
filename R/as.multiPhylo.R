#' Convert object to `multiPhylo` class
#'
#' Converts representations of phylogenetic trees to an object of the "ape"
#' class `multiPhylo`.
#' @param x Object to be converted
#' @return `as.multiPhylo` returns an object of class `multiPhylo`
#' @examples
#' as.multiPhylo(BalancedTree(8))
#' as.multiPhylo(list(BalancedTree(8), PectinateTree(8)))
#' data("Lobo")
#' as.multiPhylo(Lobo.phy)
#' @family utility functions
#' @export
as.multiPhylo <- function(x) UseMethod("as.multiPhylo")

#' @rdname as.multiPhylo
#' @export
as.multiPhylo.phylo <- function(x) c(x)


#' @rdname as.multiPhylo
#' @export
as.multiPhylo.list <- function(x) structure(x, class = "multiPhylo")

#' @rdname as.multiPhylo
#' @return `as.multiPhylo.phyDat()` returns a list of trees, each corresponding
#' to the partitions implied by each non-ambiguous character in `x`.
#' @importFrom ape read.tree
#' @export
as.multiPhylo.phyDat <- function(x) {
  at <- attributes(x)
  cont <- at[["contrast"]]
  if ("-" %fin% colnames(cont)) {
    cont[cont[, "-"] > 0, ] <- 1
  }
  ambiguous <- rowSums(cont) != 1
  labels <- names(x)

  mat <- matrix(unlist(x), length(x), byrow = TRUE)
  mat[ambiguous[mat]] <- NA
  mat <- apply(mat, 2, function(x) {
    uniques <- table(x) == 1
    x[x %fin% names(uniques[uniques])] <- NA
    x
  })
  structure(apply(mat, 2, function(split) {
    a <- !is.na(split)
    aSplit <- split[a]
    tokens <- unique(aSplit)
    aLabels <- labels[a]
    if (length(tokens) == 1L) {
      read.tree(text = paste0("(", paste(aLabels, collapse = ", "), ");"))
    } else {
      read.tree(text = paste0("((",
                              paste(vapply(unique(aSplit), function(token) {
                                paste(aLabels[aSplit == token], collapse = ", ")
                              }, character(1)), collapse = "), ("),
                              "));"))
    }
  })[at[["index"]]], class = "multiPhylo")
}

#' @rdname as.multiPhylo
#' @importFrom ape read.tree
#' @export
as.multiPhylo.Splits <- function(x) {
  labels <- TipLabels(x)
  structure(apply(as.logical(x), 1, function(a) {
    read.tree(text = paste0(
      "((",
      paste0(labels[a], collapse = ","),
      "),(",
      paste0(labels[!a], collapse = ","),
      "));"))
  }), class = "multiPhylo")
}
