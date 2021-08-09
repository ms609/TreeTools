#' Visualise position of rogue taxa
#'
#' Plots a consensus of trees with a rogue taxon omitted, with edges coloured
#' according to the proportion of trees in which the taxon attaches to that
#' edge, after Klopfstein &amp; Spasojevic (2019).
#'
#' Note: this function is currently under development and is not fully tested.
#' Please check that results are in line with expectations before using the
#' output, and [report any errors](https://github.com/ms609/TreeTools/issues/53).
#'
#' @param trees List or `multiPhylo` object containing phylogenetic trees
#' of class `phylo` to be summarized.
#' @param tip Numeric or character identifying rogue leaf, in format accepted
#' by `drop.tip()`.
#' @param p A numeric value between 0.5 and 1 giving the proportion for a clade
#' to be represented in the consensus tree (see `consensus()`).
#' @param \dots Additional parameters to `plot.phylo()`.
#' @param plot Logical specifying whether to plot the tree.
#' @param Palette Function that takes a parameter `n` and generates a colour
#' palette with `n` entries.
#' @param nullCol Colour to paint regions of the tree on which the rogue is
#' never found.
#' @param thin,fat Numeric specifying width to plot edges if the rogue tip
#' never / sometimes does attach to them.
#' @return `RoguePlot()` returns a list whose elements are:
#' - `cons`: The reduced consensus tree, in preorder;
#' - `onEdge`: a vector of integers specifying the number of
#' trees in `trees` in which the rogue leaf is attached to each edge in turn
#' of the consensus tree;
#' - `atNode`: a vector of integers specifying the number of trees in `trees`
#' in which the rogue leaf is attached to an edge collapsed into each node
#' of the consensus tree.
#' @references
#' \insertRef{Klopfstein2019}{TreeTools}
#' @examples
#' trees <- list(read.tree(text = '(a, (b, (c, (rogue, (d, (e, f))))));'),
#'               read.tree(text = '(a, (b, (c, (rogue, (d, (e, f))))));'),
#'               read.tree(text = '(a, (b, (c, (rogue, (d, (e, f))))));'),
#'               read.tree(text = '(a, (b, (c, (rogue, (d, (e, f))))));'),
#'               read.tree(text = '(rogue, (a, (b, (c, (d, (e, f))))));'),
#'               read.tree(text = '((rogue, a), (b, (c, (d, (e, f)))));'),
#'               read.tree(text = '(a, (b, ((c, d), (rogue, (e, f)))));'),
#'               read.tree(text = '(a, (b, ((c, (rogue, d)), (e, f))));'),
#'               read.tree(text = '(a, (b, (c, (d, (rogue, (e, f))))));'))
#' RoguePlot(trees, 'rogue')
#' @template MRS
#' @importFrom ape consensus
#' @importFrom fastmatch fmatch %fin%
#' @importFrom grDevices colorRampPalette
#' @importFrom phangorn allDescendants
#' @export
RoguePlot <- function (trees, tip, p = 1, plot = TRUE,
                       Palette = colorRampPalette(c(par('fg'), 'red'),
                                                  space = 'Lab'),
                       nullCol = 'lightgrey',
                       thin = par('lwd'), fat = thin + 1L,
                       ...) {
  trees <- RenumberTips(trees, trees[[1]])
  if (is.character(tip)) {
    tip <- fmatch(tip, trees[[1]]$tip.label)
  }

  noRogue <- trees
  noRogue[] <- lapply(trees, DropTip, tip)
  noRogue[] <- lapply(noRogue, AddTip, 0, 'Dummy root')
  cons <- RootTree(consensus(noRogue, p = p), 'Dummy root') # RootTree gives Preorder
  consTip <- NTip(cons)

  nTip <- NTip(trees[[1]])
  #allTips <- logical(ceiling((nTip) / 8) * 8L + 2L) # Multiple of 8 for packBits
  allTips <- logical(nTip + 1L) # including dummy
  # leaf 1 is defined as outgroup and thus always below rogue.
  aboveRogue <- vapply(trees, function (tr) {
    edge <- AddTip(tr, 0, 'Dummy Root')$edge
    parent <- edge[, 1]
    child <- edge[, 2]
    rogueEdge <- fmatch(tip, child)
    rogueParent <- parent[rogueEdge]
    aboveRogue <- fmatch(rogueParent, child)
    splitTips <- allTips
    if (is.na(aboveRogue)) {
      splitKids <- child
    } else {
      edgeInSplit <- DescendantEdges(edge = aboveRogue,
                                     parent = parent, child = child)
      splitKids <- child[edgeInSplit]
    }
    splitTips[splitKids[splitKids <= nTip + 1L]] <- TRUE
    splitTips <- splitTips[-tip]

    #ret <- packBits(splitTips[-1]) # TODO can do something clever like this?
    if (splitTips[1]) {
      !splitTips[-1]
    } else {
      splitTips[-1]
    }
  #}, raw((length(allTips) - 2L)  / 8L))
  }, logical((length(allTips) - 2L)))
  if (is.null(dim(aboveRogue))) {
    aboveRogue <- matrix(aboveRogue, nrow = 1L)
  }

  # Vector, each entry corresponds to a tree
  tipsAboveRogue <- colSums(aboveRogue)
  nAtTip <- c(sum(tipsAboveRogue == nTip - 1L), # At first tip
              double(nTip - 1L))
  atTip <- tipsAboveRogue == 1L
  tipMatches <- 1L + apply(aboveRogue[, atTip, drop = FALSE], 2, which)
  tab <- table(as.integer(tipMatches))
  nAtTip[as.integer(names(tab))] <- tab

  unmatchedTrees <- !(tipsAboveRogue %fin% c(0L, 1L, nTip - 1L))
  consSplits <- PolarizeSplits(as.Splits(cons), 1)
  splits <- !as.logical(consSplits)
  nSplits <- nrow(splits)
  edgeMatches <- match(data.frame(aboveRogue[, unmatchedTrees, drop = FALSE]),
                       data.frame(t(splits[, -1])))

  nAtSplit <- double(nSplits)
  tab <- table(edgeMatches)
  nAtSplit[as.integer(names(tab))] <- tab
  decipher <- c(nAtTip, rep.int(NA, cons$Nnode))
  decipher[as.integer(names(consSplits))] <- nAtSplit
  nOnEdge <- decipher[cons$edge[, 2]]


  nAtNode <- c(nAtTip[length(nAtTip)], double(cons$Nnode - 1L))
  unmatchedTrees[unmatchedTrees] <- is.na(edgeMatches)
  if (any(unmatchedTrees)) {
    nodeDescs <- vapply(allDescendants(cons)[-seq_len(nTip)], function (tips) {
      ret <- logical(nTip)
      ret[tips[tips <= nTip]] <- TRUE
      if (ret[1]) !ret[-1] else ret[-1]
    }, logical(nTip - 1L))
    aboveRogue <- rbind(FALSE, aboveRogue) # tip A is below rogue by definition

    nNode <- ncol(nodeDescs)
    consEdge <- cons$edge
    parentOf <- c(NA_integer_, vapply(nTip + seq_len(nNode)[-1], function (node) {
      consEdge[consEdge[, 2] == node, 1]
    }, integer(1))) - nTip


    for (i in seq_len(nNode)[-1]) {
      overlap <- aboveRogue[nodeDescs[, i], unmatchedTrees, drop = FALSE]
      active <- apply(overlap, 2, any)
      within <- apply(overlap, 2, all)
      unmatchedTrees[unmatchedTrees][active & within] <- FALSE
      nAtNode[parentOf[i]] <- nAtNode[parentOf[i]] + sum(active & within)
      if (!any(unmatchedTrees)) {
        break
      }
    }
  }

  cons <- DropTip(cons, 'Dummy root')
  nOnEdge <- nOnEdge[2:(length(nOnEdge) - 1L)]
  nAtNode <- nAtNode[-1]

  if (plot) {
    pal <- c(NA, Palette(length(trees)))
    plot.phylo(cons,
               edge.color = ifelse(nOnEdge > 0, pal[nOnEdge + 1L], nullCol),
               node.color = ifelse(nAtNode > 0,
                                   pal[c(double(consTip - 1L), nAtNode) + 1L],
                                   nullCol),
               edge.width = ifelse(nOnEdge > 0, fat, thin),
               node.width = ifelse(c(double(consTip - 1L), nAtNode) > 0,
                                   fat, thin),
               ...)
  }

  # Return:
  list(cons = cons, onEdge = nOnEdge, atNode = nAtNode)
}
