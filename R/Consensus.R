#' Construct consensus trees
#'
#' `Consensus()` calculates the consensus of a set of trees, using the
#' algorithm of \insertCite{Day1985}{TreeTools}.
#'
#' @param trees List of trees, optionally of class `multiPhylo`.
#' @param p Proportion of trees that must contain a split for it to be reported
#' in the consensus.  `p = 0.5` gives the majority-rule consensus; `p = 1` (the
#' default) gives the strict consensus.
#' @param check.labels Logical specifying whether to check that all trees have
#' identical labels.  Defaults to `TRUE`, which is slower.
#'
#' @return `Consensus()` returns an object of class `phylo`, rooted as in the
#' first entry of `trees`.
#' @examples
#' Consensus(as.phylo(0:2, 8))
#' @seealso
#' `TreeDist::ConsensusInfo()` calculates the information content of a consensus
#' tree.
#' @template MRS
#' @family consensus tree functions
#' @family tree characterization functions
#' @references
#' \insertAllCited{}
#' @export
Consensus <- function(trees, p = 1, check.labels = TRUE) {
  if (length(trees) == 1L) {
    return(trees[[1]])
  }
  if (inherits(trees, "phylo")) {
    return(trees)
  }
  if (!is.list(trees) || is.data.frame(trees)) {
    stop("Expecting `trees` to be a list.")
  }
  
  # Remove irrelevant metadata so we don't waste time processing it
  trees <- lapply(c(trees), function(tr) {
    tr[["edge.length"]] <- NULL
    tr[["node.label"]] <- NULL
    tr
  })
  
  repeat {
    nTip <- NTip(trees)
    if (length(unique(nTip)) > 1) {
      warning("Tree sizes differ; removing leaves not in smallest.")
      trees <- lapply(trees, KeepTip, trees[[which.min(nTip)]][["tip.label"]])
    } else {
      nTip <- nTip[[1]]
      break
    }
  }
  if (nTip < 4L) {
    return(trees[[1]])
  }
  if (check.labels) {
    trees <- RenumberTips(trees, trees[[1]])
  }
  if (p < 0.5 || p > 1) {
    stop("`p` must be between 0.5 and 1.")
  }
  trees <- Preorder(trees) # Per #168; could be dispensed with with further
                           # investigation of consensus_tree
  tree1 <- trees[[1]] # Must be in Preorder for DescendantEdges()
  edg <- tree1[["edge"]]
  root <- edg[DescendantEdges(edg[, 1], edg[, 2], edge = 1), 2]
  root <- root[root <= NTip(tree1)]

  # Return:
  RootTree(.PreorderTree(
    edge = splits_to_edge(consensus_tree(trees, p), nTip),
    tip.label = TipLabels(trees[[1]])
  ), root)
}

#' Reduced consensus, omitting specified taxa
#'
#' `ConsensusWithout()` displays a consensus plot with specified taxa excluded,
#' which can be a useful way to increase the resolution of a consensus tree
#' when a few wildcard taxa obscure a consistent set of relationships.
#' `MarkMissing()` adds missing taxa as loose leaves on the plot.
#'
#' @param trees A list of phylogenetic trees, of class `multiPhylo` or `list`.
#' @param tip A character vector specifying the names (or numbers) of tips to
#' drop (using `ape::drop.tip()`).
#'
#' @return `ConsensusWithout()` returns a consensus tree (of class `phylo`)
#' without the excluded taxa.
#'
#' @examples
#' oldPar <- par(mfrow = c(1, 2), mar = rep(0.5, 4))
#'
#' # Two trees differing only in placement of tip 2:
#' trees <- as.phylo(c(0, 53), 6)
#' plot(trees[[1]])
#' plot(trees[[2]])
#'
#' # Strict consensus (left panel) lacks resolution:
#' plot(ape::consensus(trees))
#'
#' # But omitting tip two (right panel) reveals shared structure in common:
#' plot(ConsensusWithout(trees, "t2"))
#' MarkMissing("t2")
#'
#' par(oldPar)
#' @family tree manipulation
#' @family tree properties
#' @family consensus tree functions
#'
#' @template MRS
#' @export
ConsensusWithout <- function(trees, tip = character(0), ...) {
  UseMethod("ConsensusWithout")
}

#' @rdname ConsensusWithout
#' @export
ConsensusWithout.phylo <- function(trees, tip = character(0), ...) {
  DropTip(trees, tip = tip)
}

#' @rdname ConsensusWithout
#' @export
ConsensusWithout.multiPhylo <- function(trees, tip = character(0), ...) {
  Consensus(lapply(trees, DropTip, tip = tip), ...)
}

#' @rdname ConsensusWithout
#' @export
ConsensusWithout.list <- ConsensusWithout.multiPhylo

#' @rdname ConsensusWithout
#' @param position Where to plot the missing taxa.
#' See [`legend()`] for options.
#' @param \dots Additional parameters to pass on to [`ape::consensus()`] or
#' [`legend()`].
#' @return `MarkMissing()` provides a null return, after plotting the specified
#' `tip`s as a legend.
#' @importFrom graphics legend
#' @export
MarkMissing <- function(tip, position = "bottomleft", ...) {                   # nocov start
  if (length(tip) > 0) {
    legend(position, legend = gsub("_", " ", tip, fixed = TRUE),
           lwd = 1, lty = 2, bty = "n", ...)
  }
}                                                                               # nocov end
