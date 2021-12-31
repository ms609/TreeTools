#' Sort tree
#'
#' `SortTree()` sorts each node into a consistent order, so that node rotation
#' does not obscure similarities between similar trees.
#'
#' At each node, clades will be listed in `tree$edge` in decreasing size order.
#'
#' Clades that contain the same number of leaves are sorted in decreasing order
#' of minimum leaf number, so (2, 3) will occur before (1, 4).
#'
#' As trees are plotted from 'bottom up', the largest clades will 'sink' to the
#' bottom of a plotted tree.
#'
#TODO:
#' `tree` must (presently) be binary ([#25](https://github.com/ms609/TreeTools/issues/25)).
#'
#' @param tree One or more trees of class `phylo`, optionally as a list
#' or a `multiPhylo` object.
#'
#' @return `SortTree()` returns tree in the format of `tree`, with each node
#' in each tree sorted such that the larger clade is first.
#'
#' @seealso `Preorder()` also rearranges trees into a consistent shape, but
#' based on the index of leaves rather than the size of subtrees.
#'
#' @examples
#' messyTree <- as.phylo(10, 6)
#' plot(messyTree)
#'
#' sorted <- SortTree(messyTree)
#' plot(sorted)
#' ape::nodelabels()
#' ape::edgelabels()
#'
#' @family tree manipulation
#'
#' @template MRS
#' @importFrom phangorn Descendants
#' @export
SortTree <- function (tree) UseMethod('SortTree')

#' @export
#' @rdname SortTree
SortTree.phylo <- function (tree) {
  edge <- tree$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  tipLabels <- tree$tip.label
  tree.ntip <- length(tipLabels)
  if (tree.ntip + tree.ntip - 2L != length(parent)) {
    stop("`tree` must be binary.")
  }

  descendants <- Descendants(tree)
  nDescendants <- vapply(descendants, length, integer(1))
  MinKid <- function (tips) min(tipLabels[tips])
  swaps <- vapply(tree.ntip + seq_len(tree$Nnode), function(node) {
    kids <- child[parent == node]
    descs <- nDescendants[kids]
    if (all(descs == 1L)) {
      order(tipLabels[kids], method = 'radix')[1] == 1
    } else if (descs[1] == descs[2]) {
      order(vapply(descendants[kids], MinKid, character(1)),
            method = 'radix')[1] == 1
    } else {
      descs[1] < descs[2]
    }
  }, logical(1))
  for (node in tree.ntip + rev(which(swaps))) {
    childEdges <- parent == node
    kids <- child[childEdges]
    child[childEdges][2:1] <- kids
  }
  tree$edge[, 1] <- parent
  tree$edge[, 2] <- child
  attr(tree, 'order') <- NULL
  Renumber(tree)
}

#' @export
#' @rdname SortTree
SortTree.list <- function (tree) lapply(tree, SortTree)

#' @export
#' @rdname SortTree
SortTree.multiPhylo <- function (tree) {
  tree[] <- SortTree.list(tree)
  tree
}
