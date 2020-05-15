#' Consensus without taxa
#'
#' Displays a consensus plot with selected taxa excluded.
#'
#' A useful way to increase the resolution of a consensus tree when a few
#' wildcard taxa obscure a consistent set of relationships.
#'
#' @param trees A list of phylogenetic trees, of class `multiPhylo` or `list`.
#' @param tip A character vector specifying the names (or numbers) of tips to
#' drop (using `ape::drop.tip`).
#' @param \dots Additional parameters to pass on to [`ape::consensus()`] or
#' [`legend()`].
#'
#' @return `ConsensusWithout` returns a consensus tree (of class `phylo`)
#' without the excluded taxa.
#'
#' @examples
#' oldPar <- par(mfrow=c(1, 2), mar=rep(0.5, 4))
#'
#' # Two trees differing only in placement of tip 2:
#' trees <- as.phylo(c(0, 53), 6)
#' plot(trees[[1]])
#' plot(trees[[2]])
#'
#' # Strict consensus lacks resolution:
#' plot(ape::consensus(trees))
#'
#' # But omitting tip two reveals shared structure in common:
#' plot(ConsensusWithout(trees, 't2'))
#' MarkMissing('t2')
#'
#' par(oldPar)
#' @family tree manipulation
#' @family tree properties
#'
#' @template MRS
#' @export
ConsensusWithout <- function (trees, tip = character(0), ...) {
  UseMethod('ConsensusWithout')
}

#' @importFrom ape consensus
#' @rdname ConsensusWithout
#' @export
ConsensusWithout.phylo <- function (trees, tip = character(0), ...) {
  DropTip(trees, tip = tip)
}

#' @importFrom ape consensus
#' @rdname ConsensusWithout
#' @export
ConsensusWithout.multiPhylo <- function (trees, tip = character(0), ...) {
  consensus(lapply(trees, DropTip, tip = tip), ...)
}

#' @rdname ConsensusWithout
#' @export
ConsensusWithout.list <- ConsensusWithout.multiPhylo

#' @describeIn ConsensusWithout Adds labels for taxa omitted from a plotted
#' consensus tree.
#' @param position Where to plot the missing taxa.  See [legend] for options.
#' @importFrom graphics legend
#' @template MRS
#' @export
MarkMissing <- function (tip, position='bottomleft', ...) {
  if (length(tip) > 0) {
    legend(position, legend=gsub('_', ' ', tip, fixed=TRUE),
         lwd=1, lty=2, bty='n', ...)
  }
}

#' Sort tree
#'
#' Sorts each node into a consistent order, so similar trees look visually
#' similar.
#'
#' @template treeParam
#'
#' @return `SortTree` returns a tree of class `phylo`, with each node sorted
#' such that the larger clade is first.
#'
#' @seealso [`RenumberTree()`]
#'
#' @family tree manipulation
#'
#' @template MRS
#' @export
SortTree <- function(tree) {
  edge <- tree$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  tipLabels <- tree$tip.label
  tree.ntip <- length(tipLabels)
  descendants <- Descendants(tree)
  nDescendants <- vapply(descendants, length, integer(1))
  MinKid <- function (tips) min(tipLabels[tips])
  swaps <- vapply(tree.ntip + seq_len(tree$Nnode), function(node) {
    kids <- child[parent == node]
    descs <- nDescendants[kids]
    if (all(descs == 1L)) {
      order(tipLabels[kids])[1] == 1
    } else if (descs[1] == descs[2]) {
      order(vapply(descendants[kids], MinKid, character(1)))[1] == 1
    } else {
      descs[1] < descs[2]
    }
  }, logical(1))
  for (node in tree.ntip + rev(which(swaps))) {
    childEdges <- parent==node
    kids <- child[childEdges]
    child[childEdges][2:1] <- kids
  }
  tree$edge[, 1] <- parent
  tree$edge[, 2] <- child
  attr(tree, 'order') <- NULL
  Cladewise(Renumber(tree))
}

#' Write Newick Tree
#'
#' Writes a tree in Newick format.  This differs from ape's `write.tree`
#' in the encoding of spaces as spaces, rather than underscores.
#'
#' @template treeParam
#'
#' @return `NewickTree` returns a character string denoting `tree` in Newick
#' format.
#'
#' @examples
#' NewickTree(BalancedTree(6))
#'
#' @seealso [`as.Newick`]
#' @importFrom ape write.tree
#' @export
NewickTree <- function(tree) gsub('_', ' ', write.tree(tree), fixed=TRUE)
