#' Remove metadata from trees
#' 
#' `TopologyOnly()` removes all information from trees except for their
#' topologies and leaf labels.  This allows other functions to process
#' trees more rapidly, as they do not need to process unneeded metadata.
#' 
#' @inheritParams Preorder
#' @return Returns `tree`, with each tree in [`Preorder`], with edge lengths,
#' node labels and other attributes removed.
#' @template MRS
#' @export
TopologyOnly <- function(tree) UseMethod("TopologyOnly")

#' @export
TopologyOnly.phylo <- function(tree) {
  Preorder(structure(list(edge = tree[["edge"]],
                          Nnode = tree[["Nnode"]],
                          tip.label = tree[["tip.label"]]),
                     order = attr(tree, "order"),
                     class = "phylo"))
}

#' @export
TopologyOnly.multiPhylo <- function(tree) {
  tree[] <- lapply(tree, TopologyOnly)
  tree
}

#' @export
TopologyOnly.list <- function(tree) {
  lapply(tree, TopologyOnly)
}

#' @export
TopologyOnly.default <- function(tree) {
  tree
}
