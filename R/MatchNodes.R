#' Match nodes and edges between trees
#' 
#' `MatchNodes()` and `MatchEdges()` matches nodes or edges in one tree to
#' entries in the second that denote a clade with identical tip labels.
#' 
#' `MatchNodes()` treats clades as identical only if they contain exactly the
#' same tips.
#' `MatchEdges()` disregards tips of `table` that are absent from `x`.
#' 
#' @param x Tree whose nodes are to be matched.
#' @param table Tree containing nodes to be matched against.
#' @param tips Logical specifying whether to return matches for tips;
#' unless `TRUE`, only the matches for internal nodes will be returned.
#' @return `MatchNodes()` returns an integer vector giving, for each internal
#' node of `x` (and, if `tips = TRUE`, each tip, listed first), the number of
#' the node in `table` that subtends the same clade, or `nomatch` if no node
#' corresponds.
#' 
#' `MatchEdges()` returns an integer vector giving, for each edge of `x`,
#' the index of the edge in `table` that subtends the same clade, or `nomatch`.
#' Where several nodes or edges of `table` denote the same clade (because of
#' nodes with a single child), the first is returned.
#' @examples
#' MatchNodes(BalancedTree(8), RootTree(BalancedTree(8)))
#' MatchEdges(BalancedTree(8), RootTree(BalancedTree(8)))
#' @inheritParams match
#' @template MRS
#' @family tree navigation
#' @family tree properties
#' @export
MatchEdges <- function(x, table, nomatch = NA_integer_) {
  xEdge <- .Edge(x)
  tableEdge <- .Edge(table)
  xLab <- TipLabels(x)
  tabLab <- TipLabels(table)
  tipMatch <- if (is.null(xLab) && is.null(tabLab)) {
    # Match tips by number
    xTip <- min(xEdge[, 1]) - 1
    tabTip <- min(tableEdge[, 1]) - 1
    `length<-`(seq_len(tabTip), xTip)
  } else if (is.null(xLab) || is.null(tabLab) ||
             anyDuplicated(xLab) || anyDuplicated(tabLab)) {
    return(.MatchEdgesR(xEdge, tableEdge, xLab, tabLab, nomatch))
  } else {
    match(xLab, tabLab)
  }
  ret <- match_clades(xEdge, tableEdge, tipMatch, TRUE)
  
  # Return:
  `[<-`(ret, is.na(ret), nomatch)
}

# Reference implementation; retained for inputs with missing or duplicated
# tip labels, whose semantics the linear-time algorithm does not reproduce.
.MatchEdgesR <- function(xEdge, tableEdge, xLab, tabLab,
                         nomatch = NA_integer_) {
  tipMatch <- if (is.null(xLab) && is.null(tabLab)) {
    # Match tips by number
    xTip <- min(xEdge[, 1]) - 1
    tabTip <- min(tableEdge[, 1]) - 1
    `length<-`(seq_len(tabTip), xTip)
  } else {
    match(xLab, tabLab)
  }
  seek <- DescendantTips(xEdge[, 1], xEdge[, 2])
  find <- DescendantTips(tableEdge[, 1], tableEdge[, 2])[, tipMatch]
  find[is.na(find)] <- FALSE
  ret <- match(as.data.frame(t(seek)), as.data.frame(t(find)))
  
  # Return:
  `[<-`(ret, is.na(ret), nomatch)
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
  if (inherits(x, "phylo")) {
    xEdge <- x[["edge"]]
    xLab <- x[["tip.label"]]
  } else {
    xEdge <- .Edge(x)
    xLab <- TipLabels(x)
  }
  if (inherits(table, "phylo")) {
    tableEdge <- table[["edge"]]
    tabLab <- table[["tip.label"]]
  } else {
    tableEdge <- .Edge(table)
    tabLab <- TipLabels(table)
  }
  if (is.null(xLab) || is.null(tabLab) ||
      anyDuplicated(xLab) || anyDuplicated(tabLab)) {
    return(.MatchNodesR(xEdge, tableEdge, xLab, tabLab, nomatch, tips))
  }
  ret <- match_clades(xEdge, tableEdge, match(xLab, tabLab), FALSE)
  if (!isTRUE(tips)) {
    ret <- ret[-seq_along(xLab)]
  }
  
  # Return:
  `[<-`(ret, is.na(ret), nomatch)
}

# Reference implementation; see .MatchEdgesR()
.MatchNodesR <- function(xEdge, tableEdge, xLab, tabLab,
                         nomatch = NA_integer_, tips = FALSE) {
  allLab <- union(xLab, tabLab)
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
                                     newTips = TipLabels(old), ...) {
  nodeLabel[MatchNodes(list(edge = new, tip.label = newTips),
             old, tips = FALSE) - NTip(old)]
}

#' @export
.UpdateNodeLabel.phylo <- function(new, old, nodeLabel = old[["node.label"]], ...) {
  nodeLabel[MatchNodes(new, old, tips = FALSE) - NTip(old)]
}
