#' Sort tree
#'
#' `SortTree()` sorts each node into a consistent order, so that node rotation
#' does not obscure similarities between similar trees.
#'
#' At each node, clades will be listed in `tree[["edge"]]` in decreasing size
#' order.
#'
#' Clades that contain the same number of leaves are sorted in decreasing order
#' of minimum leaf number, so (2, 3) will occur before (1, 4).
#'
#' As trees are plotted from 'bottom up', the largest clades will 'sink' to the
#' bottom of a plotted tree.
#'
#' @param tree One or more trees of class `phylo`, optionally as a list
#' or a `multiPhylo` object.
#' @param how Character vector specifying sort method:
#' `"Cladesize"` rotates each node such that the larger clade is first,
#' thus appearing lower when plotted;
#' `"TipLabels"` rotates nodes such that labels listed sooner in `order`
#' are listed first, and thus plot lower.
#' @param order Character vector listing tip labels in sequence they should
#' appear on tree. Clades containing a taxon earlier in this list will be listed
#' sooner and thus plot lower on a tree.  Taxa not listed in `order` will be
#' treated as if they were last in the list.
#'
#' @return `SortTree()` returns tree in the format of `tree`, with each node
#' in each tree sorted 
#'
#' @seealso `Preorder()` also rearranges trees into a consistent shape, but
#' based on the index of leaves rather than the size of subtrees.
#'
#' @examples
#' messyTree <- as.phylo(10, 6)
#' plot(messyTree)
#'
#' sorted <- SortTree(messyTree)
#' plot(sorted)
#' ape::nodelabels()
#' ape::edgelabels()
#' ape::tiplabels(adj = c(2, 1/3))
#' 
#' plot(SortTree(messyTree, how = "tip"))
#' @template MRS
#' @family tree manipulation
#'
#' @export
SortTree <- function(tree, how = "cladesize", order = TipLabels(tree)) {
  UseMethod('SortTree')
}

#' @export
#' @rdname SortTree
SortTree.phylo <- function(tree, how = "cladesize", order = TipLabels(tree)) {
  edge <- tree[["edge"]]
  parent <- edge[, 1]
  child <- edge[, 2]
  tipLabels <- tree[["tip.label"]]
  nTip <- length(tipLabels)
  if (!is.null(tree[["edge.length"]])) {
    warning("Edge lengths are not supported (#49)")
    tree[["edge.length"]] <- NULL
  }
  
  descendants <- .ListDescendents(tree)
  method <- pmatch(tolower(how), c("cladesize", "tiplabels"))
  if (is.na(method)) {
    warning("`how` must be an unambiguous contraction of ",
              "\"cladesize\" or \"tiplabels\"")
    method <- 1L
  }
    
  
  switch(method,
    { # Clade Size
      nDescendants <- vapply(descendants, length, integer(1))
      MinKid <- function(tips) min(tipLabels[tips])
      for (node in nTip + seq_len(tree[["Nnode"]])) {
        childEdges <- parent == node
        kids <- child[childEdges]
        newOrder <- order(nDescendants[kids],
                          vapply(descendants[kids], MinKid, tipLabels[1]),
                          method = "radix", decreasing = TRUE)
        child[childEdges] <- kids[newOrder]
      }
    }, { # Tip labels
      weights <- match(TipLabels(tree), order, nomatch = NTip(tree) + 1L)
      MinKid <- function(tips) min(weights[tips])
      for (node in nTip + seq_len(tree[["Nnode"]])) {
        childEdges <- parent == node
        kids <- child[childEdges]
        newOrder <- order(vapply(descendants[kids], MinKid, integer(1)),
                          method = "radix")
        child[childEdges] <- kids[newOrder]
    }
  })
  
 
  tree[["edge"]][, 2] <- child
  attr(tree, 'order') <- NULL
  Renumber(tree)
}

.ListDescendents <- function(tree) {
  edge <- Postorder(tree[["edge"]])
  parent <- edge[, 1]
  child <- edge[, 2]
  # Every node occurs once in `child` except the root
  descendants <- vector('list', length(child) + 1L)
  descendants[seq_len(NTip(tree))] <- seq_len(NTip(tree))
  for (i in seq_along(parent)) {
    descendants[[parent[i]]] <- c(descendants[[parent[i]]], 
                                  descendants[[child[i]]])
  }
  
  # Return:
  descendants
}

#' @export
#' @rdname SortTree
SortTree.list <- function(tree) lapply(tree, SortTree)

#' @export
#' @rdname SortTree
SortTree.multiPhylo <- function(tree) {
  tree[] <- SortTree.list(tree)
  tree
}
