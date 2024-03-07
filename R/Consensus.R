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
Consensus <- function(trees, p = 1, check.labels = TRUE, info = FALSE) {
  if (length(trees) == 1L) {
    return(trees[[1]])
  }
  if (inherits(trees, "phylo")) {
    return(trees)
  }
  if (!is.list(trees) || is.data.frame(trees)) {
    stop("Expecting `trees` to be a list.")
  }
  repeat {
    nTip <- NTip(trees)
    if (length(unique(nTip)) > 1) {
      warning("Tree sizes differ; removing leaves not in smallest.")
      trees <- lapply(trees, KeepTip, trees[[which.min(nTip)]][["tip.label"]])
    } else {
      nTip <- nTip[1]
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
  if (isTRUE(info)) {
    info <- "phylogenetic"
  }
  
  if (!is.logical(info)) {
    mode <- pmatch(tolower(info),
                  c("phylogenetic", "clustering", "credibility"))
    if (is.na(mode)) {
      stop("`info` must be 'phylogenetic', 'clustering' or 'credibility'")
    }
    counts <- .CountSplits(trees)
    duplicate <- duplicated(counts[["splits"]])
    allSplits <- as.Splits(counts[["splits"]][!duplicate, ],
                           tipLabels = TipLabels(trees[[1]]))
    info <- counts[[switch(mode, "pic", "cic", "count")]][!duplicate]
    kept <- logical(length(info))
    active <- !kept
    repeat {
      best <- which.max(info * active)
      kept[best] <- TRUE
      active[best] <- FALSE
      compatible <- CompatibleSplits(allSplits[[best]], allSplits[[active]])
      active[active][!compatible] <- FALSE
      if (!any(active)) {
        break
      }
    }
    splits <- allSplits[[kept]]
  } else {
    splits <- as.Splits(consensus_tree(trees, p),
                        tipLabels = TipLabels(trees[[1]]))
  }
  tree1 <- Preorder(trees[[1]])
  edg <- tree1[["edge"]]
  root <- edg[DescendantEdges(edg[, 1], edg[, 2], edge = 1), 2]
  root <- root[root <= NTip(tree1)]

  # Return:
  RootTree(as.phylo(splits), root)
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

#' Count splits and calculate information content
#' 
#' Counts the number of times that each split occurs in a posterior sample of
#' trees, and computes the information content implied by split membership
#' and probability.
#' 
#' This is marked as an internal function as its behaviour is subject to change;
#' it may be replaced with more efficient functions for specific applications
#' in future. Please alert the maintainer if you plan to call the function
#' directly in your own work.
#' 
#' @param trees List or multiPhylo object containing phylogenetic trees of 
#' class "phylo", on the same set of tips.
#' @return `.CountSplits()` currently returns a list with four elements:
#' - `splits` is a logical matrix; each column corresponds to a leaf, each row
#' to a split
#' - `count` lists the number of occurrences of that split
#' - `pic` lists the phylogenetic information content of that split, assuming
#' that split frequency is proportional to split probability
#' - `cic` lists the clustering information content of each split, on the same
#' basis.
#' Splits are listed tree by tree, including duplicates.
#' The order of splits corresponds to their position in the cluster table
#' representation of each tree, so may differ from the sequence given by
#' `as.Splits(trees)`.
#' @template MRS
#' @keywords internal
#' @export
.CountSplits <- function(trees) {
  if (!is.null(attr(trees, "TipLabel"))) {
    # Decompress tip labels
    trees <- lapply(trees, I)
  }
  if (!all(vapply(trees, inherits, TRUE, "phylo"))) {
    stop("All `trees` must be of class \"phylo\"")
  }
  if (length(unique(NTip(trees))) > 1) {
    stop("All `trees` must contain same leaves")
  }
  
  # Return:
  count_splits(trees)
}

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
