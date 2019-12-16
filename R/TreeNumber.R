#' Unique integer indices for bifurcating tree topologies
#'
#' Functions converting between phylogenetic trees and their unique decimal
#' representation.
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
#' 0 &times; (1 &times; 3 &times; 5 &times; 9) +
#'
#' 0 &times; (1 &times; 3 &times; 5) +
#'
#' 3 &times; (1 &times; 3) +
#'
#' 1 &times; (1)
#'
#' = 10
#'
#' Note that the hyperexponential nature of tree space means that there are &gt;
#' 2^30 unique 12-tip trees.  As integers &gt; 2^31 are not supported by R,
#' numbers representing larger trees are represented internally as a vector of
#' nine-digit integer 'chunks'.
#'
#' @param x Integer identifying the tree (see details).
#' @param nTip Integer specifying number of tips in the tree.
#' @template tipLabelsParam
#' @param \dots Additional parameters for consistency with S3 methods (unused).
#'
#' @return `as.phylo.numeric` returns a tree of class `phylo`.
#'
#' @examples
#' tree <- as.phylo(10, nTip = 6)
#' plot(tree)
#' as.TreeNumber(tree)
#'
#' # Larger trees:
#' as.TreeNumber(BalancedTree(19))
#'
#' # If > 9 digits, represent the tree number as a string.
#' treeNumber <- as.TreeNumber("1234567890123", nTip = 14)
#' tree <- as.phylo(treeNumber)
#'
#' @exportClass TreeNumber
#' @name TreeNumber
#

#' @rdname TreeNumber
#'
#' @return `as.TreeNumber` returns an object of class `TreeNumber`,
#' which comprises a numeric vector, whose elements represent successive
#' nine-digit chunks of the decimal integer corresponding to the tree topology
#' (in big endian order).  The `TreeNumber` object has attributes
#' `nTip` and `tip.labels`.
#' @export
as.TreeNumber <- function(x, ...) UseMethod('as.TreeNumber')

#' @rdname TreeNumber
#' @importFrom ape root
#' @export
as.TreeNumber.phylo <- function (x, ...) {
  x <- root(x, 1, resolve.root = TRUE)
  edge <- x$edge
  nTip <- NTip(x)
  edge <- PostorderEdges(edge[, 1], edge[, 2], nTip = nTip)
  structure(edge_to_num(edge[[1]], edge[[2]], nTip),
            nTip = nTip,
            tip.labels = TipLabels(x),
            class = 'TreeNumber')
}

#' @rdname TreeNumber
#' @export
as.TreeNumber.character <- function (x, nTip, tipLabels = TipLabels(nTip), ...) {
  len <- nchar(x)
  ends <- rev(len - (seq_len((len / 9L) + 1L) - 1L) * 9L)
  starts <- pmax(1L, ends - 8L)

  # Return:
  structure(as.integer(substring(x, starts, ends)),
            nTip = nTip,
            tip.labels = tipLabels,
            class = 'TreeNumber')
}

#' @rdname TreeNumber
#' @template MRS
#' @references Based on a concept by John Tromp, employed in Li _et al._ 1996.
#'
#' \insertRef{Li1996}{TreeTools}
#'
#' @importFrom ape as.phylo
#' @export
as.phylo.numeric <- function (x, nTip = attr(x, 'nTip'),
                              tipLabels = attr(x, 'tip.label'), ...) {
  if (is.null(tipLabels)) tipLabels <- paste0('t', seq_len(nTip))
  edge <- RenumberEdges(num_to_parent(x, nTip), seq_len(nTip + nTip - 2L))
  structure(list(edge = do.call(cbind, edge),
                 tip.label = tipLabels,
                 Nnode = nTip - 1L),
            order = 'postorder',
            class = 'phylo')
}

#' @rdname TreeNumber
#' @export
as.phylo.TreeNumber <- function (x, ...) as.phylo.numeric(x, ...)

#' Print TreeNumber object
#'
#' S3 method for objects of class `TreeNumber`.
#'
#' @param x Object of class `TreeNumber`.
#' @param \dots Additional arguments for consistency with S3 method (unused).
#'
#' @export
print.TreeNumber <- function (x, ...) {
  nTip <- attr(x, 'nTip')
  cat("Phylogenetic tree number", .PrintedTreeNumber(x), "of", NUnrooted(nTip),
      "\n", nTip, "tips:", paste0(attr(x, 'tip.labels')))
}

#' @keywords internal
.PrintedTreeNumber <- function (x) {
  paste0(c(x[1], formatC(x[-1], digits = 9L, flag='0')), collapse = '')
}
