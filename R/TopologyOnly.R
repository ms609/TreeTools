#' Remove metadata from trees
#' 
#' `TopologyOnly()` removes all information from trees except for their
#' topologies and leaf labels.
#' 
#' @inheritParams Preorder
#' @return Returns `tree`, with each tree in [`Preorder`], with edge lengths,
#' node labels and other attributes removed.
#' @template MRS
#' @export
TopologyOnly <- function(tree) UseMethod("TopologyOnly")

#' @export
TopologyOnly.phylo <- function(tree) {
  startOrder <- attr(tree, "order")
  ret <- structure(list(edge = tree[["edge"]],
                        Nnode = tree[["Nnode"]],
                        tip.label = tree[["tip.label"]]),
                   order = startOrder,
                   class = "phylo")
  if (length(startOrder) && startOrder == "preorder") {
    ret
  } else {
    Preorder(ret)
  }
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
