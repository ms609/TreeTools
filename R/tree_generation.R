#' Extract tip names from a dataset of unknown class
#'
#' @template tipsForTreeGeneration
#'
#' @return Character vector listing tip names.
#'
#' @template MRS
#' @keywords internal
GetTipNames <- function (tips) {
  if (mode(tips) == 'numeric') {
    if (length(tips) == 1L) {
      tips <- paste0('t', seq_len(tips))
    } else {
      tips <- as.character(tips)
    }
  } else if (class(tips) == 'phyDat') {
    tips <- names(tips)
  } else if (class(tips) == 'phylo') {
    tips <- tips$tip.label
  }

  # Return:
  tips
}

#' Generate random tree topology
#'
#' @template tipsForTreeGeneration
#' @param root Taxon to use as root (if desired; FALSE otherwise)
#'
#' @template MRS
#' @importFrom ape rtree
#' @importFrom ape root
#' @family tree generation functions
#' @export
RandomTree <- function (tips, root = FALSE) {
  tips <- GetTipNames(tips)
  nTips <- length(tips)
  tree <- rtree(nTips, tip.label=tips, br=NULL)
  if (root != FALSE) {
    tree <- root(tree, root, resolve.root=TRUE)
  }

  # Return:
  tree
}

#' Generate a Pectinate Tree
#'
#' Generates a pectinate (caterpillar) tree with the specified tip labels.
#'
#' @template tipsForTreeGeneration
#'
#' @return A tree of class `phylo`.
#' @family tree generation functions
#' @template MRS
#' @export
PectinateTree <- function (tips) {
  tips <- GetTipNames(tips)
  nTips <- length(tips)

  nEdge <- nTips + nTips - 2L
  tipSeq <- seq_len(nTips - 1L)

  parent <- rep(seq_len(nTips - 1L) + nTips, each = 2L)

  child <- integer(nEdge)
  child[tipSeq + tipSeq - 1L] <- tipSeq
  child[tipSeq + tipSeq] <- tipSeq + nTips + 1L
  child[nEdge] <- nTips

  structure(list(
    edge = matrix(c(parent, child), ncol = 2L),
    Nnode = nTips - 1L,
    tip.label = tips
  ), order = 'cladewise', class='phylo')
}

#' Generate a Balanced Tree
#'
#' Generates a balanced (symmetrical) tree with the specified tip labels.
#'
#' @template tipsForTreeGeneration
#'
#' @return A tree of class `phylo`.
#' @family tree generation functions
#' @template MRS
#' @importFrom ape read.tree
#' @export
BalancedTree <- function (tips) {
  tips <- GetTipNames(tips)

  # Return:
  read.tree(text=paste0(BalancedBit(tips), ';'))
}

BalancedBit <- function (tips, nTips = length(tips)) {
  if (nTips < 4L) {
    if (nTips == 2L) {
      paste0('(', tips[1L], ',', tips[2L], ')')
    } else if (nTips == 3L) {
      paste0('(', tips[1L], ',(', tips[2L], ',', tips[3L], '))')
    } else {
      tips
    }
  } else {
    # Recurse:
    firstHalf <- seq_len(nTips / 2L)
    paste0('(', BalancedBit(tips[firstHalf]), ',',
           BalancedBit(tips[seq_along(tips)[-firstHalf]]), ')')
  }
}

#' Neighbour Joining Tree
#'
#' Generates a rooted neighbour joining tree, with no edge lengths
#'
#' @template datasetParam
#'
#' @return an object of class \code{phylo}
#'
#' @template MRS
#' @importFrom ape nj root
#' @importFrom phangorn dist.hamming
#' @family tree generation functions
#' @export
NJTree <- function (dataset) {
  nj.tree <- nj(dist.hamming(dataset))
  nj.tree <- root(nj.tree, outgroup=names(dataset)[1], resolve.root=TRUE)
  nj.tree$edge.length <- NULL
  nj.tree
}

#' Force taxa to form an outgroup
#'
#' Given a tree or a list of taxa, rearrange the ingroup and outgroup taxa such that the two
#' are sister taxa across the root, without altering the relationships within the ingroup
#' or within the outgroup.
#'
#' @param tree either a tree of class \code{phylo}, or a character vector listing the names of
#'        all the taxa in the tree, from which a random tree will be generated.
#' @param outgroup a vector containing the names of taxa to include in the outgroup
#'
#' @return a tree where all outgroup taxa are sister to all remaining taxa,
#'         otherwise retaining the topology of the ingroup.
#' @template MRS
#' @importFrom ape rtree
#' @importFrom ape root drop.tip bind.tree
#' @export
EnforceOutgroup <- function (tree, outgroup) {
  if (class(tree) == 'phylo') {
    taxa <- tree$tip.label
  } else if (class(tree) == 'character') {
    tree <- root(rtree(length(taxa), tip.label=taxa, br=NULL), taxa[1], resolve.root=TRUE)
  } else {
    stop ("tree must be of class phylo")
  }

  if (length(outgroup) == 1) return (root(tree, outgroup, resolve.root=TRUE))

  ingroup <- taxa[!(taxa %in% outgroup)]
  if (!all(outgroup %in% taxa) || length(ingroup) + length(outgroup) != length(taxa)) {
    stop ("All outgroup taxa must occur in speficied taxa")
  }

  ingroup.branch <- drop.tip(tree, outgroup)
  outgroup.branch <- drop.tip(tree, ingroup)

  result <- root(bind.tree(outgroup.branch, ingroup.branch, 0, 1), outgroup, resolve.root=TRUE)
  RenumberTips(Renumber(result), taxa)
}

#' Unique integer indices for bifurcating tree topologies
#'
#' There are `NUnrooted(n)` unrooted trees with _n_ tips.
#' As such, each _n_-tip tree can be uniquely identifyed by a non-negative
#' integer _x_ < `NUnrooted(n)`.
#'
#' This integer can be converted by a tree by treating it as a mixed-base number,
#' with bases 1, 3, 5, 7, ... (2_n_ - 5).
#'
#' Each digit of this mixed base number corresponds to a tip, and determines
#' the location on a growing tree to which that tip should be added.
#'
#' We start with a two-tip tree, and treat 0 as the origin of the tree.
#'
#' ```
#'
#' 0 ---- 1
#'
#' ```
#'
#' We add tip 2 by breaking an edge and inserting a node (numbered 2 + nTip - 1).
#' In this example, we'll work up to a six-tip tree; this node will be numbered
#' 2 + 6 - 1 = 7.
#' There is only one edge on which tip 2 can be added.  Let's add node 7 and tip 2:
#'
#' ```
#'
#' 0 ---- 7 ---- 1
#'        |
#'        |
#'        2
#'
#' ```
#'
#' There are now three edges on which tip 3 can be added.  Our options are:
#' Option 0: the edge leading to 1;
#' Option 1: the edge leading to 2;
#' Option 2: the edge leading to 7.
#'
#' If we select option 1, we produce:
#'
#' ```
#'
#' 0 ---- 7 ---- 1
#'        |
#'        |
#'        8 ---- 2
#'        |
#'        |
#'        3
#'
#' ```
#' `1` is now the final digit of our mixed-base number
#'
#' There are five places to add tip 4:
#' Option 0: the edge leading to 1;
#' Option 1: the edge leading to 2;
#' Option 2: the edge leading to 3;
#' Option 3: the edge leading to 7;
#' Option 4: the edge leading to 8.
#'
#' If we chose option 3, then `3` would be the penultimate digit of our
#' mixed-base number
#'
#' If we chose option 0 for the next two additions, we could specify this tree
#' with the mixed-base number 0021.  We can convert this into decimal:
#'
#' 0 × (1 × 3 × 5 × 9) +
#' 0 × (1 × 3 × 5) +
#' 3 × (1 × 3) +
#' 1 × (1)
#' = 10
#'
#' Note that the hyperexponential nature of tree space means that there are
#' > 2^30 unique 12-tip trees.  As
#'
#' @param x Integer identifying the tree (see details).
#' @param nTip Integer specifying number of tips in the tree.
#' @template tipLabelsParam
#' @template MRS
#' @references Based on a concept by John Tromp (1995)
#' @importFrom ape as.phylo
#' @export
as.phylo.numeric <- function (x, nTip = attr(x, 'nTip'),
                              tipLabels = attr(x, 'tip.label')) {
  if (is.null(tipLabels)) tipLabels <- paste0('t', seq_len(nTip))
  edge <- RenumberEdges(num_to_parent(x, nTip), seq_len(nTip + nTip - 2L))
  structure(list(edge = do.call(cbind, edge),
                 tip.label = tipLabels,
                 Nnode = nTip - 1L),
            order = 'postorder',
            class = 'phylo')
}

#' @describeIn as.phylo.numeric Converts tree to index.
#' @importFrom ape root
#' @export
as.numeric.phylo <- function (x) {
  x <- root(x, 1, resolve.root = TRUE)
  edge <- x$edge
  nTip <- NTip(x)
  edge <- PostorderEdges(edge[, 1], edge[, 2], nTip = nTip)
  structure(edge_to_num(edge[[1]], edge[[2]], nTip),
            nTip = nTip,
            tip.labels = TipLabels(x))
}
