#' Visualize position of rogue taxa
#'
#' Plots a consensus of trees with a rogue taxon omitted, with edges coloured
#' according to the proportion of trees in which the taxon attaches to that
#' edge, after \insertCite{Klopfstein2019;textual}{TreeTools}.
#' 
#' Rogue taxa can be identified using the package \pkg{Rogue}
#' \insertCite{SmithCons}{TreeTools}.
#'
#' @param trees List or `multiPhylo` object containing phylogenetic trees
#' of class `phylo` to be summarized.
#' @param tip Numeric or character identifying rogue leaf, in format accepted
#' by [`DropTip()`].
#' @param p A numeric value between 0.5 and 1 giving the proportion for a clade
#' to be represented in the consensus tree (see [`Consensus()`]).
#' @param \dots Additional parameters to `plot.phylo()`.
#' @param plot Logical specifying whether to plot the tree.
#' @param Palette Function that takes a parameter `n` and generates a colour
#' palette with `n` entries.
#' @param nullCol Colour to paint regions of the tree on which the rogue is
#' never found.
#' @template edgeLengthParam
#' @inheritParams RootTree
#' @param thin,fat Numeric specifying width to plot edges if the rogue tip
#' never / sometimes does attach to them.
#' @param sort Logical specifying whether to sort consensus tree using
#' [`SortTree()`].
#' @param legend Character vector specifying position of legend (e.g.
#' `"bottomleft"`), or `"none"` to suppress legend.
#' For fine-grained control of legend, use [`PlotTools::SpectrumLegend()`].
#' @param legend.inset Numeric specifying fraction of plot width / height
#' by which the legend's position should be inset.
#'
#' @return `RoguePlot()` invisibly returns a list whose elements are:
#' - `cons`: The reduced consensus tree, in preorder;
#' - `onEdge`: a vector of integers specifying the number of
#' trees in `trees` in which the rogue leaf is attached to each edge in turn
#' of the consensus tree;
#' - `atNode`: a vector of integers specifying the number of trees in `trees`
#' in which the rogue leaf is attached to an edge collapsed into each node
#' of the consensus tree.
#' - `legendLabels`: A character vector suggesting labels for a plot legend;
#' suitable for `PlotTools::SpectrumLegend(legend = x$legendLabels)`.
#' @references \insertAllCited{}
#' @examples
#' trees <- list(read.tree(text = "(a, (b, (c, (rogue, (d, (e, f))))));"),
#'               read.tree(text = "(a, (b, (c, (rogue, (d, (e, f))))));"),
#'               read.tree(text = "(a, (b, (c, (rogue, (d, (e, f))))));"),
#'               read.tree(text = "(a, (b, (c, (rogue, (d, (e, f))))));"),
#'               read.tree(text = "(rogue, (a, (b, (c, (d, (e, f))))));"),
#'               read.tree(text = "((rogue, a), (b, (c, (d, (e, f)))));"),
#'               read.tree(text = "(a, (b, ((c, d), (rogue, (e, f)))));"),
#'               read.tree(text = "(a, (b, ((c, (rogue, d)), (e, f))));"),
#'               read.tree(text = "(a, (b, (c, (d, (rogue, (e, f))))));"))
#' plotted <- RoguePlot(trees, "rogue", legend = "topleft", legend.inset = 0.02)
#' PlotTools::SpectrumLegend(
#'   "bottomleft",
#'   palette = colorRampPalette(c(par("fg"), "#009E73"), space = "Lab")(100),
#'   legend = plotted$legendLabels,
#'   cex = 0.4
#' )
#' @template MRS
#' @importFrom fastmatch fmatch %fin%
#' @importFrom graphics par
#' @importFrom grDevices colorRamp colorRampPalette rgb
#' @importFrom PlotTools SpectrumLegend
#' @family consensus tree functions
#' @export
RoguePlot <- function(trees, tip, p = 1, plot = TRUE,
                      Palette = colorRampPalette(c(par("fg"), "#009E73"),
                                                 space = "Lab"),
                      nullCol = rgb(colorRamp(unlist(par(c("fg", "bg"))),
                                              space = "Lab")(0.8) / 255),
                      edgeLength = NULL,
                      thin = par("lwd"), fat = thin + 1L,
                      outgroupTips,
                      sort = FALSE,
                      legend = "none",
                      legend.inset = 0,
                      ...) {
  tipLabels <- TipLabels(trees[[1]])
  nTip <- length(tipLabels)

  at <- attributes(trees)
  if (!missing(outgroupTips)) {
    trees <- lapply(trees, RootTree, intersect(outgroupTips, tipLabels))
  }
  trees <- RenumberTips(trees, tipLabels)
  attributes(trees) <- at

  noRogue <- trees
  attr(noRogue, "TipLabel") <- NULL
  dummyRoot <- "xxTREETOOLSxxDUMMYxxROOTxx"
  noRogue[] <- lapply(lapply(noRogue, DropTip, tip),
                      AddTip, 0, dummyRoot)
  class(noRogue) <- "multiPhylo"
  cons <- RootTree(Consensus(noRogue, p = p, check.labels = FALSE),
                   dummyRoot) # RootTree gives Preorder
  if (sort) {
    cons <- SortTree(cons)
  }
  consTip <- NTip(cons)

  if (is.character(tip)) {
    tip <- fmatch(tip, tipLabels)
  }
  pole <- nTip

  #allTips <- logical(ceiling((nTip) / 8) * 8L + 2L) # Multiple of 8 for packBits
  allTips <- logical(nTip + 1L) # including dummy
  # dummyRoot is, by definition, always below rogue.
  aboveRogue <- .vapply(trees, function(tr) {
    edge <- AddTip(tr, 0, dummyRoot)[["edge"]]
    parent <- edge[, 1]
    child <- edge[, 2]
    rogueEdge <- fmatch(tip, child)
    rogueParent <- parent[[rogueEdge]]
    aboveRogue <- fmatch(rogueParent, child) # edge above parent of rogue tip
    splitTips <- allTips
    stopifnot(!is.na(aboveRogue))
    edgeInSplit <- DescendantEdges(parent = parent, child = child,
                                   edge = aboveRogue)
    splitKids <- child[edgeInSplit]
    splitTips[splitKids[splitKids <= nTip + 1L]] <- TRUE
    splitTips <- splitTips[-tip]

    #ret <- packBits(splitTips[-1]) # TODO can do something clever like this?
    
    # `pole` is the dummy root, so cannot have been in the rogue's sister clade
    # We thus guarantee the polarity of splitTips
    stopifnot(!splitTips[[pole]]) # if NOT, return !splitTips[-pole]
    splitTips[-pole]
  #}, raw((length(allTips) - 2L)  / 8L))
  }, logical((length(allTips) - 2L)))

  # Vector, each entry corresponds to a tree
  tipsAboveRogue <- colSums(aboveRogue) # includes dummy root
  
  # Initialize count of trees where rogue occurs alongside a leaf
  nAtTip <- c(double(nTip - 1L),
              sum(tipsAboveRogue == nTip - 1L)) # At pole
  # Populate nAtTip
  atTip <- tipsAboveRogue == 1L
  tipMatches <- apply(aboveRogue[, atTip, drop = FALSE], 2, which)
  tab <- tabulate(as.integer(tipMatches))
  nAtTip[seq_along(tab)] <- tab

  unmatchedTrees <- !(tipsAboveRogue %fin% c(0L, 1L, nTip - 1L))
  consSplits <- PolarizeSplits(as.Splits(
    cons, tipLabels = c(tipLabels[-tip], dummyRoot)), pole)
  splits <- !as.logical(consSplits)
  nSplits <- nrow(splits)
  # Had previously used fmatch here, but encountered unexpected error
  # Error in fmatch: cannot take a writable DATAPTR of an ALTLIST
  #  [class: wrap_list, pkg: base]
  # Replaced with match for safety (?)
  edgeMatches <- match(data.frame(aboveRogue[, unmatchedTrees, drop = FALSE]),
                       data.frame(t(splits[, -pole, drop = FALSE])))

  nAtSplit <- double(nSplits)
  tab <- tabulate(edgeMatches)
  nAtSplit[which(as.logical(tab))] <- tab[as.logical(tab)]
  decipher <- c(nAtTip, rep.int(NA, cons[["Nnode"]]))
  decipher[as.integer(names(consSplits))] <- nAtSplit
  nOnEdge <- decipher[cons[["edge"]][, 2]]


  nAtNode <- c(nAtTip[length(nAtTip)], double(cons[["Nnode"]] - 2L))
  unmatchedTrees[unmatchedTrees] <- is.na(edgeMatches)
  if (any(unmatchedTrees)) {
    nodeDescs <- .NodeDescendants(cons, nTip, pole)
    nNode <- nrow(nodeDescs)

    unmatchedGroup <- aboveRogue[, unmatchedTrees, drop = FALSE]
    nUnmatched <- sum(unmatchedTrees)

    # nodeOverlapsGroup <- apply(unmatchedGroup, 2, function(gp) {
    #   apply(nodeDescs[gp, , drop = FALSE], 2, any)
    # })
    floors <- .apply(unmatchedGroup, 2, function(gp) {
       nodeContainsAllGroup <- apply(nodeDescs[, gp, drop = FALSE], 1, all)
       nodeContainsNonGroup <- apply(nodeDescs[, !gp, drop = FALSE], 1, any)
       nodeContainsAllGroup & nodeContainsNonGroup
    })

    # active[i, ] != within[i, ], or we'd be on an edge
    atNode <- table(apply(cbind(floors), 2,
                          function(x) 1 + length(x) - which.max(rev(x))))
    nAtNode[as.integer(names(atNode))] <- atNode

  }

  cons <- DropTip(cons, dummyRoot)
  if (!is.null(edgeLength)) {
    cons[["edge.length"]] <- rep_len(edgeLength, dim(cons[["edge"]])[1])
  }
  nOnEdge <- nOnEdge[2:(length(nOnEdge) - 1L)]
  maxVal <- max(c(nOnEdge, nAtNode)) + 1L

  legendLabels <- if (maxVal < 8) {
    maxVal:1
  } else if ((maxVal - 1) %% 5 == 0) {
    seq(maxVal, 1, length.out = 6)
  } else if ((maxVal - 1) %% 4 == 0) {
    seq(maxVal, 1, length.out = 5)
  } else if ((maxVal - 1) %% 3 == 0) {
    seq(maxVal, 1, length.out = 4)
  } else if ((maxVal - 1) %% 6 == 0) {
    seq(maxVal, 1, length.out = 7)
  } else {
    ceiling(seq(maxVal, 1, length.out = 5))
  }
  legendLabels <- c(paste(legendLabels[-length(legendLabels)], "trees"),
                    "1 tree")
  
  if (plot) {
    #pal <- c(NA, Palette(length(trees)))
    pal <- Palette(maxVal)
    plot(cons,
         edge.color = ifelse(nOnEdge > 0, pal[nOnEdge + 1L], nullCol),
         node.color = c(double(consTip - 1L),
                        ifelse(nAtNode > 0, pal[nAtNode + 1L], nullCol)),
         edge.width = ifelse(nOnEdge > 0, fat, thin),
         node.width = ifelse(c(double(consTip - 1L), nAtNode) > 0,
                             fat, thin),
         ...)
    if (legend != "none") {
      PlotTools::SpectrumLegend(
        legend,
        bty = "n",
        palette = pal,
        legend = legendLabels,
        inset = legend.inset
      )
    }
  }

  # Return:
  invisible(list(cons = cons, onEdge = nOnEdge, atNode = nAtNode,
                 legendLabels = legendLabels))
}

# `tree` must be in preorder
.NodeDescendants <- function(tree, nTip, pole) {
  edge <- tree[["edge"]]
  parent <- edge[, 1]
  child <- edge[, 2]
  descs <- matrix(FALSE, max(parent), nTip)
  diag(descs[seq_len(nTip), ]) <- TRUE
  
  for (i in rev(seq_along(parent))) {
    descs[parent[i], ] <- descs[parent[i], ] | descs[child[i], ]
  }
  
  # Internal nodes only
  descs <- descs[-seq_len(nTip), , drop = FALSE]
  
  # Ignore the dummy root node
  descs <- descs[-1, , drop = FALSE]
  
  flip <- descs[, pole]
  descs[flip, ] <- !descs[flip, ]
  
  # Return:
  descs[, -pole, drop = FALSE]
}

# TODO in r4.1.0 use apply(simplify = false)?
.vapply <- function(...) {
  x <- vapply(...)
  if (is.null(dim(x))) {
    x <- matrix(x, 1L)
  }

  # Return:
  x
}

.apply <- function(...) {
  x <- apply(...)
  if (is.null(dim(x))) {
    x <- matrix(x, 1L)
  }

  # Return:
  x
}
