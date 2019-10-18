#' Generate random tree topology from dataset
#' 
#' @param dataset A dataset in \code{\link[phangorn]{phyDat}} format
#' @param root Taxon to use as root (if desired; FALSE otherwise)
#' 
#' @author Martin R. Smith 
#' @importFrom ape rtree
#' @importFrom ape root
#' @family tree generation functions
#' @export
RandomTree <- function (dataset, root = FALSE) {
  tree <- rtree(length(dataset), tip.label=names(dataset), br=NULL)
  return (if (root != FALSE) root(tree, root, resolve.root=TRUE) else tree)
}

#' Neighbour Joining Tree
#' 
#' Generates a rooted neighbour joining tree, with no edge lengths
#'
#' @template datasetParam
#' 
#' @return an object of class \code{phylo}
#'
#' @author Martin R. Smith
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
#' @author Martin R. Smith
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

#' Generate a Pectinate Tree
#' 
#' Generates a pectinate (caterpillar) tree with the specified tips
#' 
#' @param tips An integer specifying the number of tips, or a character vector
#' naming the tips, or a dataset phylogenetic dataset whose names will be used
#' to name the tips.
#' 
#' @return A tree of class `phylo`.
#' @family tree generation functions
#' @author Martin R. Smith
#' @export
PectinateTree <- function (tips) {
  if (length(tips) == 1L && mode(tips) == 'numeric') {
    nTips <- tips
    tips <- paste0('t', seq_len(nTips))
  } else if (class(tips) == 'phyDat') {
    nTips <- length(tips)
    tips <- names(tips)
  } else {
    nTips <- length(tips)
  }
  
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
