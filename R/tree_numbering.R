#' Reorder edges of a phylogenetic tree
#'
#' Wrappers for the C functions called by
#' \code{ape::\link[ape:reorder.phylo]{reorder.phylo}}.
#' These call the C functions directly, so are faster -- but don't perform
#' as many checks on user input.  Bad input could crash R.
#'
#'
#' @param nTip,nNode,nEdge Integer specifying the number of tips, nodes
#' and edges in the input tree.
#' @template treeParent
#' @template treeChild
#'
#' @return `NeworderPruningwise` returns an integer vector specifying the
#' pruningwise order of edges within a tree.
#'
#' @examples
#' nTip <- 8L
#' tree <- BalancedTree(nTip)
#' edge <- tree[["edge"]]
#' pruningwise <- NeworderPruningwise(nTip, tree$Nnode, edge[, 1], edge[, 2],
#'                                    dim(edge)[1])
#' cladewise <- NeworderPhylo(nTip, edge[, 1], edge[, 2], dim(edge)[1], 1L)
#' postorder <- NeworderPhylo(nTip, edge[, 1], edge[, 2], dim(edge)[1], 2L)
#'
#' tree[["edge"]] <- tree[["edge"]][pruningwise, ]
#'
#' @author
#'  - C algorithm: Emmanuel Paradis
#'  - R wrapper: Martin R. Smith
#' @family C wrappers
#' @keywords internal
#' @useDynLib TreeTools, .registration = TRUE
#' @name Neworder
#' @export
NeworderPruningwise <- function(nTip, nNode, parent, child, nEdge) {
  .C('ape_neworder_pruningwise', as.integer(nTip), as.integer(nNode),
     as.integer(parent), as.integer(child), as.integer(nEdge),
     integer(nEdge))[[6]]
}

#' @rdname Neworder
#' @param `whichwise` Integer specifying whether to order edges (1)
#' cladewise; or (2) in postorder.
#' @return `NeworderPhylo` returns an integer vector specifying the order
#' of edges under the ordering sequence specified by `whichwise`.
#' @keywords internal
#' @export
NeworderPhylo <- function(nTip, parent, child, nEdge, whichwise) {
  .C('ape_neworder_phylo', as.integer(nTip), as.integer(parent),
     as.integer(child), as.integer(nEdge), integer(nEdge),
     as.integer(whichwise), NAOK = TRUE)[[5]]
}

#' @rdname Reorder
#'
#' @template treeParent
#' @template treeChild
#' @param weight Optional vector specifying the weight of each edge,
#' corresponding to the `edge.length` property of a `phylo` object.
#'
#' @return `RenumberTree()` returns an edge matrix for a tree of class `phylo`
#' following the preorder convention for edge and node numbering.
#'
#' @family tree manipulation
#' @family C wrappers
#' @export
RenumberTree <- function(parent, child, weight) {
  if (missing(weight)) {
    .Call(`_TreeTools_preorder_edges_and_nodes`, parent, child)
  } else {
    .Call(`_TreeTools_preorder_weighted`, parent, child, weight)
  }
}

#' @rdname Reorder
#'
#' @param \dots Deprecated; included for compatibility with previous versions.
#' @return `RenumberEdges()` formats the output of `RenumberTree()` into a list
#' whose two entries correspond to the new parent and child vectors.
#' @export
RenumberEdges <- function(parent, child, ...) {
  oenn <- .Call(`_TreeTools_preorder_edges_and_nodes`, parent, child)

  # Return:
  list(oenn[, 1], oenn[, 2])
}

#' Reorder tree edges and nodes
#'
#' Functions for systematically ordering the internal edges of trees.
#'
#' `Reorder()` is a wrapper for \code{ape:::.reorder_ape}.
#' Calling this C function directly is approximately twice as fast as using
#' \code{ape::\link[ape:reorder.phylo]{cladewise}} or
#' \code{ape::\link[ape:reorder.phylo]{postorder}}
#'
#' `Cladewise()`, `ApePostorder()` and `Pruningwise()` are convenience
#' functions to the corresponding functions in 'ape'. Single nodes may
#' need to be collapsed using [ape::collapse.singles] first.  'ape' functions
#' can cause crashes if nodes are numbered unconventionally -- sometimes
#' encountered after using tree rearrangement functions, e.g. `phangorn::SPR()`.
#'
#' `Preorder()` is more robust: it supports polytomies, nodes may be numbered
#' in any sequence, and edges may be listed in any order in the input tree.
#' Its output is guaranteed to be identical for any tree of an equivalent
#' topology, allowing unique trees to be detected by comparing sorted edge
#' matrices alone.
#'
#' Nodes and edges in a preorder tree are numbered starting from the deepest
#' node.
#' Each node is numbered in the sequence in which it is encountered, and
#' each edge is listed in the sequence in which it is visited.
#'
#' At each node, child edges are sorted from left to right in order of the
#' lowest-numbered leaf in the subtree subtended by each edge; i.e. an edge
#' that leads eventually to tip 1 will be to the left of an edge leading to a
#' subtree containing tip 2.
#'
#' Numbering begins by following the leftmost edge of the root node,
#' and sorting its descendant subtree into preorder.
#' Then, the next edge at the root node is followed, and its descendants
#' sorted into preorder, until each edge has been visited.
#'
#' `RenumberTree()` and `RenumberEdges()` are wrappers for the C function
#' `preorder_edges_and_nodes()`; they do not perform the same checks on input
#' as `Preorder()` and are intended for use where performance is at a premium.
#'
#'
#' `Postorder()` numbers nodes as in `Preorder()`, and lists edges in
#' descending order of parent node number.  If a tree is already
#' in postorder, it will not be rearranged unless `force = TRUE`.
#'
#' Methods applied to numeric inputs do not check input for sanity, so should
#' be used with caution: malformed input may cause undefined results, including
#' crashing R.
#'
#' Trees with >8191 leaves require additional memory and are not handled
#' by `Postorder()` at present.
#' If you need to process such large trees, please contact the maintainer for
#' advice.
#'
#' @template treeParam
#' @template nTipParam
#' @param edge Two-column matrix listing the parent and child of each edge in a
#' tree, corresponding to `tree[["edge"]]`. Optional in `Cladewise()`.
#' @param renumber Logical specifying whether to renumber nodes such that they
#' increase in number away from the root.
#'
#' @return `ApePostorder()`, `Cladewise()`, `Postorder()`, `Preorder()` and
#' `Pruningwise()` each return a tree of class `phylo` with nodes following the
#' specified numbering scheme.
#' @author
#' `Preorder()` and `Postorder()`: Martin R. Smith.
#'
#' `Cladewise()`, `ApePostorder()` and `Pruningwise()`: modified by Martin R.
#' Smith from \code{.reorder_ape()} in \pkg{ape} (Emmanuel Paradis).
#'
#'
#' @seealso Rotate each node into a consistent orientation with [`SortTree()`].
#'
#' @family C wrappers
#' @keywords internal
#' @name Reorder
NULL

#' @describeIn Reorder Reorder tree cladewise.
#' @export
Cladewise <- function(tree, nTip, edge) UseMethod('Cladewise')

#' @rdname Reorder
#' @export
Cladewise.phylo <- function(tree, nTip = NTip(tree), edge = tree[["edge"]]) {
  if (!is.null(attr(tree, "order")) && attr(tree, "order") == "cladewise") {
    return(tree)
  }
  nEdge <- dim(edge)[1]
  nNode <- tree[["Nnode"]]
  if (nNode == 1) {
    return(tree)
  }
  if (nNode >= nTip) {
    stop("`tree` apparently badly conformed")
  }

  newOrder <- NeworderPhylo(nTip, edge[, 1], edge[, 2], nEdge, 1)

  tree[["edge"]] <- edge[newOrder, ]
  if (!is.null(tree[["edge.length"]])) {
    tree[["edge.length"]] <- tree[["edge.length"]][newOrder]
  }
  attr(tree, "order") <- "cladewise"
  attr(tree, "suborder") <- NULL

  # Return:
  tree
}

#' @rdname Reorder
#' @export
Cladewise.list <- function(tree, nTip, edge) {
  lapply(tree, Cladewise)
}

#' @rdname Reorder
#' @export
Cladewise.multiPhylo <- function(tree, nTip, edge) {
  tree[] <- lapply(tree, Cladewise)
  attr(tree, 'order') <- 'cladewise'
  attr(tree, "suborder") <- NULL
  tree
}

#' @rdname Reorder
#' @export
Cladewise.matrix <- function(tree, nTip = min(tree[, 1]) - 1L, edge) {
  if (is.numeric(tree)) {
    newOrder <- NeworderPhylo(nTip, tree[, 1], tree[, 2], dim(tree)[1], 1L)

    # Return:
    tree[newOrder, ]
  } else {
    NextMethod()
  }
}

#' @rdname Reorder
#' @export
Cladewise.NULL <- function(tree, nTip = min(tree[, 1]) - 1L, edge) NULL

#' @describeIn Reorder Reorder tree in Postorder using ape's `postorder`
#' function, which is robust to unconventional node numbering.
#' @export
ApePostorder <- function(tree, nTip, edge) UseMethod('ApePostorder')

#' @rdname Reorder
#' @export
ApePostorder.phylo <- function(tree, nTip = NTip(tree), edge = tree[["edge"]]) {
  if (!is.null(attr(tree, "order")) && attr(tree, "order") == "postorder") {
    return(tree)
  }
  nEdge <- dim(edge)[1]
  nNode <- tree[["Nnode"]]
  if (nNode == 1) {
    return(tree)
  }
  if (nNode >= nTip) {
    stop("`tree` apparently badly conformed")
  }
  neworder <- NeworderPhylo(nTip, edge[, 1], edge[, 2], nEdge, 2L)
  tree[["edge"]] <- edge[neworder, ]
  if (!is.null(tree[["edge.length"]])) {
    tree[["edge.length"]] <- tree[["edge.length"]][neworder]
  }
  attr(tree, "order") <- "postorder"
  attr(tree, "suborder") <- "ape"
  tree
}

#' @rdname Reorder
#' @export
ApePostorder.list <- function(tree, nTip, edge) {
  lapply(tree, ApePostorder)
}

#' @rdname Reorder
#' @export
ApePostorder.NULL <- function(tree, nTip, edge) NULL

#' @rdname Reorder
#' @export
ApePostorder.multiPhylo <- function(tree, nTip, edge) {
  tree[] <- lapply(tree, ApePostorder)
  attr(tree, 'order') <- 'postorder'
  attr(tree, "suborder") <- "ape"
  tree
}

#' @rdname Reorder
#' @param force Logical specifying whether to rearrange trees already in
#' postorder, in order to ensure edges are ordered in the 'TreeTools' fashion.
#' @export
Postorder <- function(tree, force = FALSE) {
  UseMethod('Postorder')
}

#' @rdname Reorder
#' @export
Postorder.phylo <- function(tree, force = FALSE) {
  if (is.null(attr(tree, "order"))
      || attr(tree, "order") != "postorder"
      || (force &&
          (is.null(attr(tree, 'suborder')) ||
           attr(tree, 'suborder') != 'TreeTools'))) {
    weight <- tree[["edge.length"]]
    edge <- tree[["edge"]]
    parent <- edge[, 1]
    child <- edge[, 2]
    if (!is.null(weight)) {
      rnt <- RenumberTree(parent, child, weight)
      edge <- rnt[[1]]
      newOrder <- order(edge[, 1], edge[, 2], decreasing = TRUE,
                        method = "radix")
      tree[["edge.length"]] <- rnt[[2]][newOrder]
    } else {
      edge <- RenumberTree(parent, child)
      newOrder <- order(edge[, 1], edge[, 2], decreasing = TRUE,
                        method = "radix")
    }
    tree[["edge"]] <- edge[newOrder, , drop = FALSE]
    attr(tree, "order") <- "postorder"
    attr(tree, "suborder") <- "TreeTools"
  }
  tree
}

#' @rdname Reorder
#' @export
Postorder.NULL <- function(tree, force = FALSE) NULL

#' @rdname Reorder
#' @export
Postorder.list <- function(tree, force = FALSE) {
  lapply(tree, Postorder, force = force)
}


#' @rdname Reorder
#' @export
Postorder.multiPhylo <- function(tree, force = FALSE) {
  tree[] <- lapply(tree, Postorder, force = force)
  attr(tree, 'order') <- 'postorder'
  tree
}

#' @rdname Reorder
#' @return `Postorder.numeric` accepts a numeric matrix corresponding to the
#' `edge` entry of a tree of class `phylo`, and returns a two-column array
#' corresponding to `tree`, with edges listed in postorder
#' @export
Postorder.numeric <- function(tree, force = FALSE) {
  edge <- RenumberTree(tree[, 1], tree[, 2])
  ordr <- order(edge[, 1], edge[, 2], decreasing = TRUE, method = "radix")
  edge[ordr, , drop = FALSE]
}

#' @describeIn Reorder Reorder tree Pruningwise.
#' @export
Pruningwise <- function(tree, nTip, edge) UseMethod('Pruningwise')

#' @rdname Reorder
#' @export
Pruningwise.phylo <- function(tree, nTip = NTip(tree),
                              edge = tree[["edge"]]) {
  if (!is.null(attr(tree, "order")) && attr(tree, "order") == 'pruningwise') {
    return(tree)
  }
  nEdge <- dim(edge)[1]
  nNode <- tree[["Nnode"]]
  if (nNode == 1) {
    return(tree)
  }
  if (nNode >= nTip) {
    stop("`tree` apparently badly conformed")
  }
  tree <- Cladewise(tree, nTip, edge)
  edge <- tree[["edge"]]
  neworder <- NeworderPruningwise(nTip, nNode, edge[, 1], edge[, 2], nEdge)
  tree[["edge"]] <- tree[["edge"]][neworder, ]
  if (!is.null(tree[["edge.length"]])) {
    tree[["edge.length"]] <- tree[["edge.length"]][neworder]
  }
  attr(tree, "order") <- "pruningwise"
  attr(tree, "suborder") <- NULL
  tree
}

#' @rdname Reorder
#' @export
Pruningwise.list <- function(tree, nTip, edge) {
  lapply(tree, Pruningwise)
}

#' @rdname Reorder
#' @export
Pruningwise.multiPhylo <- function(tree, nTip, edge) {
  tree[] <- lapply(tree, Pruningwise)
  attr(tree, 'order') <- 'pruningwise'
  attr(tree, "suborder") <- NULL
  tree
}

#' @rdname Reorder
#' @export
Pruningwise.NULL <- function(tree, nTip, edge) NULL


#' @describeIn Reorder Reorder tree in Preorder (special case of cladewise).
#' @export
Preorder <- function(tree) UseMethod('Preorder')

#' @rdname Reorder
#' @export
Preorder.phylo <- function(tree) {
  startOrder <- attr(tree, 'order')
  if (length(startOrder) && startOrder == 'preorder') {
    tree
  } else {
    edge <- tree[["edge"]]
    parent <- edge[, 1]
    child <- edge[, 2]
    lengths <- tree[["edge.length"]]
    if (is.null(lengths)) {
      tree[["edge"]] <- RenumberTree(parent, child)
    } else {
      edge <- RenumberTree(parent, child, lengths)
      tree[["edge"]] <- edge[[1]]
      tree[["edge.length"]] <- edge[[2]]
    }
    attr(tree, 'order') <- 'preorder'
    attr(tree, "suborder") <- NULL

    # Return:
    tree
  }
}

#' @rdname Reorder
#' @export
Preorder.numeric <- function(tree) {
  RenumberTree(tree[, 1], tree[, 2])
}

#' @rdname Reorder
#' @export
Preorder.multiPhylo <- function(tree) {
  tree[] <- lapply(tree, Preorder)
  attr(tree, 'order') <- 'preorder'
  attr(tree, "suborder") <- NULL
  tree
}

#' @rdname Reorder
#' @export
Preorder.list <- function(tree) {
  lapply(tree, Preorder)
}

#' @rdname Reorder
#' @export
Preorder.NULL <- function(tree) NULL


#' Renumber a tree's tips
#'
#' `RenumberTips(tree, tipOrder)` sorts the tips of a phylogenetic tree `tree`
#' such that the indices in `tree[["edge"]][, 2]` correspond to the order of
#' leaves given in `tipOrder`.
#'
#' @template treeParam
#' @param tipOrder A character vector containing the values of
#'        `tree[["tip.label"]]` in the desired sort order, or an object
#'        (perhaps of class `phylo` or `Splits`) with tip labels.
#'
#' @return `RenumberTips()` returns `tree`, with the tips' internal
#' representation numbered to match `tipOrder`.
#'
#' @examples
#' data('Lobo') # Loads the phyDat object Lobo.phy
#' tree <- RandomTree(Lobo.phy)
#' tree <- RenumberTips(tree, names(Lobo.phy))
#'
#' @family tree manipulation
#'
#' @template MRS
#' @export
RenumberTips <- function(tree, tipOrder) UseMethod('RenumberTips')

#' @rdname RenumberTips
#' @export
RenumberTips.phylo <- function(tree, tipOrder) {
  startOrder <- tree[["tip.label"]]
  newOrder <- TipLabels(tipOrder, single = TRUE)
  if (identical(startOrder, newOrder)) return(tree)
  if (length(startOrder) != length(newOrder)) {
    startOnly <- setdiff(startOrder, newOrder)
    newOnly <- setdiff(newOrder, startOrder)
    stop("Tree labels and tipOrder must match.",
         if (length(newOnly)) "\n  Missing in `tree`: ",
         paste0(newOnly, collapse = ', '),
         if (length(startOnly)) "\n  Missing in `tipOrder`: ",
         paste0(startOnly, collapse = ', ')
         )
  }

  nTip <- length(startOrder)
  child <- tree[["edge"]][, 2]
  tips <- child <= nTip

  matchOrder <- match(startOrder, newOrder)
  if (any(is.na(matchOrder))) {
    stop("Tree labels ", paste0(startOrder[is.na(matchOrder)], collapse = ', '),
         " missing from `tipOrder`")
  }
  tree[["edge"]][tips, 2] <- matchOrder[tree[["edge"]][tips, 2]]
  tree[["tip.label"]] <- newOrder
  tree
}

#' @rdname RenumberTips
#' @export
RenumberTips.multiPhylo <- function(tree, tipOrder) {
  at <- attributes(tree)
  labelled <- !is.null(at[['TipLabel']])
  tree <- lapply(tree, RenumberTips.phylo, tipOrder)
  if (labelled) {
    at[["TipLabel"]] <- TipLabels(tipOrder)
  }
  attributes(tree) <- at
  tree
}

#' @rdname RenumberTips
#' @export
RenumberTips.list <- function(tree, tipOrder) {
  lapply(tree, RenumberTips, tipOrder)
}

#' @rdname RenumberTips
#' @export
RenumberTips.NULL <- function(tree, tipOrder) NULL
