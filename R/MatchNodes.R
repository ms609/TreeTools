#' Match nodes and edges between trees
#' 
#' `MatchNodes()` and `MatchEdges()` matches nodes or edges in one tree to
#' entries in the second that denote a clade with identical tip labels.
#' 
#' The current implementation is potentially inefficient.
#' Please contact the maintainer to request a more efficient implementation if
#' this function is proving a bottleneck.
#' 
#' @param x Tree whose nodes are to be matched.
#' @param table Tree containing nodes to be matched against.
#' @examples
#' MatchNodes(BalancedTree(8), RootTree(BalancedTree(8)))
#' @inheritParams match
#' @template MRS
#' @family tree navigation
#' @family tree properties
#' @export
MatchEdges <- function(x, table, nomatch = NA_integer_) {
  xEdge <- x[["edge"]]
  tableEdge <- table[["edge"]]
  tipMatch <- match(TipLabels(x), TipLabels(table))
  seek <- DescendantTips(xEdge[, 1], xEdge[, 2])
  find <- DescendantTips(tableEdge[, 1], tableEdge[, 2])[, tipMatch]
  find[is.na(find)] <- FALSE
  match(as.data.frame(t(seek)), as.data.frame(t(find)))
}

MatchNodes <- function(x, table, nomatch = NA_integer_) {
  xEdge <- x[["edge"]]
  xLab <- TipLabels(x)
  tableEdge <- table[["edge"]]
  tabLab <- TipLabels(table)
  allLab <- union(xLab, tabLab)
  tipMatch <- match(xLab, tabLab)
  seek <- DescendantTips(xEdge[, 1], xEdge[, 2])[, match(allLab, xLab)]
  seek[is.na(seek)] <- FALSE
  find <- DescendantTips(tableEdge[, 1], tableEdge[, 2])[, match(allLab, tabLab)]
  find[is.na(find)] <- FALSE
  findFrame <-as.data.frame(cbind(t(find), allLab %in% tabLab),
                            optional = TRUE)
  matching <- match(as.data.frame(t(seek), optional = TRUE), findFrame)
  matchRoot <- match(as.data.frame(allLab %in% xLab, optional = TRUE),
                     findFrame)
  
  xRoot <- NTip(x) + 1L
  ret <- tableEdge[matching, 2][order(c(xEdge[, 2], xRoot))]
  ret[xRoot] <- tableEdge[matchRoot, 2]
  ret
}

