#' Generate random tree topology
#'
#' @template tipsForTreeGeneration
#' @param root Taxon to use as root (if desired; FALSE otherwise)
#'
#' @return `RandomTree` returns a random tree of class `phylo`, with the
#'  specified tips, and no branch lengths specified.
#'
#' @template MRS
#' @importFrom ape rtree
#' @importFrom ape root
#' @family tree generation functions
#' @export
RandomTree <- function (tips, root = FALSE) {
  tips <- TipLabels(tips)
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
  tips <- TipLabels(tips)
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
  tips <- TipLabels(tips)

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
