#' Generate a random tree topology
#'
#' Generates a binary tree with a random topology on specified tips, optionally
#' rooting the tree on a given tip.
#'
#' @template tipsForTreeGeneration
#' @param root Tip to use as root (if desired; FALSE otherwise)
#'
#' @return `RandomTree` returns a random tree of class `phylo`, with the
#'  specified tips, and no branch lengths specified.
#'
#' @examples
#' RandomTree(letters[1:5])
#'
#' data('Lobo')
#' RandomTree(Lobo.phy)
#'
#' @template MRS
#' @importFrom ape rtree root
#' @family tree generation functions
#' @export
RandomTree <- function (tips, root = FALSE) {
  tips <- TipLabels(tips)
  nTips <- length(tips)
  tree <- rtree(nTips, rooted = root != FALSE, tip.label = tips, br = NULL)
  if (root != FALSE) {
    tree <- root(tree, root, resolve.root = TRUE)
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
#' @return `PectinateTree` and `BalancedTree` each return a binary tree of
#'  class `phylo` of the specified shape.
#' @family tree generation functions
#' @examples
#' plot(PectinateTree(LETTERS[1:10]))
#'
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
  ), order = 'preorder', class = 'phylo')
}

#' Generate a Balanced Tree
#'
#' Generates a balanced (symmetrical) binary tree with the specified tip labels.
#'
#' @template tipsForTreeGeneration
#'
#' @return A tree of class `phylo`.
#' @family tree generation functions
#' @template MRS
#' @importFrom ape read.tree
#'
#' @examples
#' plot(BalancedTree(LETTERS[1:10]))

#' @export
BalancedTree <- function (tips) {
  tips <- TipLabels(tips)

  # Return:
  read.tree(text = paste0(BalancedBit(tips), ';'))
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

#' Generate a neighbour joining tree
#'
#' Generates a rooted neighbour joining tree, with no edge lengths.
#'
#' @template datasetParam
#'
#' @return `NJTree` returns an object of class \code{phylo}.
#'
#' @examples
#' data('Lobo')
#' NJTree(Lobo.phy)
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
#' Given a tree or a list of taxa, rearrange the ingroup and outgroup taxa such
#' that the two are sister taxa across the root, without changing the
#' relationships within the ingroup or within the outgroup.
#'
#' @param tree Either: a tree of class \code{phylo}; or a character vector
#' listing the names of all the taxa in the tree, from which a random tree will
#' be generated.
#' @param outgroup Character vector containing the names of taxa to include in the
#' outgroup.
#'
#' @return `EnforceOutgroup` returns a tree of class `phylo` where all outgroup
#' taxa are sister to all remaining taxa, without modifying the ingroup
#' topology.
#'
#' @examples
#' tree <- EnforceOutgroup(letters[1:9], letters[1:3])
#' plot(tree)
#'
#' @template MRS
#' @family tree manipulation
#' @export
EnforceOutgroup <- function (tree, outgroup) UseMethod('EnforceOutgroup')

#' @importFrom ape root drop.tip bind.tree
.EnforceOutgroup <- function (tree, outgroup, taxa) {
  if (length(outgroup) == 1) return (root(tree, outgroup, resolve.root = TRUE))

  ingroup <- taxa[!(taxa %in% outgroup)]
  if (!all(outgroup %in% taxa) || length(ingroup) + length(outgroup) != length(taxa)) {
    stop ("All outgroup taxa must occur in speficied taxa")
  }

  ingroup.branch <- drop.tip(tree, outgroup)
  outgroup.branch <- drop.tip(tree, ingroup)

  result <- root(bind.tree(outgroup.branch, ingroup.branch, 0, 1),
                 outgroup, resolve.root = TRUE)
  RenumberTips(Renumber(result), taxa)
}

#' @rdname EnforceOutgroup
#' @export
EnforceOutgroup.phylo <- function (tree, outgroup) {
  .EnforceOutgroup(tree, outgroup, tree$tip.label)
}

#' @rdname EnforceOutgroup
#' @importFrom ape root rtree
#' @export
EnforceOutgroup.character <- function (tree, outgroup) {
  taxa <- tree
  .EnforceOutgroup(root(rtree(length(taxa), tip.label = taxa, br = NULL),
                        taxa[1], resolve.root = TRUE),
                   outgroup, taxa)
}
