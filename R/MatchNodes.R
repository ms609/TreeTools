#' Match nodes and edges between trees
#' 
#' `MatchNodes()` and `MatchEdges()` matches nodes or edges in one tree to
#' entries in the second that denote a clade with identical tip labels.
#' 
#' Clades are compared using the cluster table of
#' \insertCite{Day1985;textual}{TreeTools}, which takes time linear in the
#' number of tips.
#' Tips are matched by label if both trees have tip labels, and otherwise by
#' number.  Tip labels must be unique.
#' 
#' `MatchNodes()` treats clades as identical only if they contain exactly the
#' same tips, so a node of `x` that subtends a tip absent from `table` has no
#' match.
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
#' `MatchEdges()` returns an integer vector giving, for each edge of `x`,
#' the index of the edge in `table` that subtends the same clade, or `nomatch`.
#' Where several nodes or edges of `table` denote the same clade (because of
#' nodes with a single child), the first is returned.
#' @examples
#' MatchNodes(BalancedTree(8), RootTree(BalancedTree(8)))
#' @references \insertAllCited{}
#' @inheritParams match
#' @template MRS
#' @family tree navigation
#' @family tree properties
#' @export
MatchEdges <- function(x, table, nomatch = NA_integer_) {
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
  ret <- match_clades(xEdge, tableEdge,
                      .TipMatch(xLab, tabLab, xEdge, tableEdge), TRUE)
  
  # Return:
  `[<-`(ret, is.na(ret), nomatch)
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
  tipMatch <- .TipMatch(xLab, tabLab, xEdge, tableEdge)
  ret <- match_clades(xEdge, tableEdge, tipMatch, FALSE)
  if (!isTRUE(tips)) {
    ret <- ret[-seq_along(tipMatch)]
  }
  
  # Return:
  `[<-`(ret, is.na(ret), nomatch)
}

# For each tip of `x`, the corresponding tip of `table`: by label if both
# trees have tip labels, otherwise by number.
.TipMatch <- function(xLab, tabLab, xEdge, tableEdge) {
  if (is.null(xLab) || is.null(tabLab)) {
    `length<-`(seq_len(.NTipOf(tabLab, tableEdge)), .NTipOf(xLab, xEdge))
  } else {
    if (anyDuplicated(xLab) || anyDuplicated(tabLab)) {
      stop("Tip labels must be unique")
    }
    match(xLab, tabLab)
  }
}

.NTipOf <- function(labels, edge) {
  if (!is.null(labels)) {
    length(labels)
  } else if (length(edge)) {
    min(edge[, 1]) - 1L
  } else {
    1L
  }
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

.UpdateNodeLabel <- function(new, old, nodeLabel = old[["node.label"]], ...) {
  UseMethod(".UpdateNodeLabel")
}

#' @export
.UpdateNodeLabel.numeric <- function(new, old, nodeLabel = old[["node.label"]],
                                     newTips = NULL, ...) {
  if (is.null(newTips)) {
    # `new` shares the tip numbering of `old`
    nodeLabel[node_label_order(new, old[["edge"]], NTip(old))]
  } else {
    nodeLabel[MatchNodes(list(edge = new, tip.label = newTips),
                         old, tips = FALSE) - NTip(old)]
  }
}

#' @export
.UpdateNodeLabel.phylo <- function(new, old, nodeLabel = old[["node.label"]], ...) {
  nodeLabel[MatchNodes(new, old, tips = FALSE) - NTip(old)]
}
