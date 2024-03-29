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
#' @param tips Logical specifying whether to return matches for tips;
#' unless `TRUE`, only the matches for internal nodes will be returned.
#' @examples
#' MatchNodes(BalancedTree(8), RootTree(BalancedTree(8)))
#' @inheritParams match
#' @template MRS
#' @family tree navigation
#' @family tree properties
#' @export
MatchEdges <- function(x, table, nomatch = NA_integer_) {
  xEdge <- .Edge(x)
  tableEdge <- .Edge(table)
  tipMatch <- match(TipLabels(x), TipLabels(table))
  seek <- DescendantTips(xEdge[, 1], xEdge[, 2])
  find <- DescendantTips(tableEdge[, 1], tableEdge[, 2])[, tipMatch]
  find[is.na(find)] <- FALSE
  match(as.data.frame(t(seek)), as.data.frame(t(find)))
}

.Edge <- function(x) {
  UseMethod(".Edge")
}

#' @export
.Edge.phylo <- function(x) {
  x[["edge"]]
}

#' @export
.Edge.numeric <- function(x) {
  x
}

#' @export
.Edge.list <- function(x) {
  if (is.null(x[["edge"]])) {
    cbind(x[[1]], x[[2]])
  } else {
    x[["edge"]]
  }
}

#' @rdname MatchEdges
#' @export
MatchNodes <- function(x, table, nomatch = NA_integer_, tips = FALSE) {
  xEdge <- .Edge(x)
  xLab <- TipLabels(x)
  tableEdge <- .Edge(table)
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
  
  nodeIndex <- c(tableEdge[, 2], length(tabLab) + 1)
  xRoot <- length(xLab) + 1L
  ret <- nodeIndex[matching][order(c(xEdge[, 2], xRoot))]
  ret[xRoot] <- nodeIndex[matchRoot]
  if (!isTRUE(tips)) {
    ret <- ret[-seq_along(xLab)]
  }
  
  # Return:
  `[<-`(ret, is.na(ret), nomatch)
}

.UpdateNodeLabel <- function(new, old, nodeLabel = old[["node.label"]], ...) {
  UseMethod(".UpdateNodeLabel")
}

#' @export
.UpdateNodeLabel.numeric <- function(new, old, nodeLabel = old[["node.label"]],
                                     newTips = TipLabels(old)) {
  nodeLabel[MatchNodes(list(edge = new, tip.label = newTips),
             old, tips = FALSE) - NTip(old)]
}

#' @export
.UpdateNodeLabel.phylo <- function(new, old, nodeLabel = old[["node.label"]]) {
  nodeLabel[MatchNodes(new, old, tips = FALSE) - NTip(old)]
}
