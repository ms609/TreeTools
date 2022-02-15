#' Unique integer indices for bifurcating tree topologies
#'
#' Functions converting between phylogenetic trees and their unique decimal
#' representation, based on a concept by John Tromp, employed in 
#' \insertCite{Li1996}{TreeTools}.
#'
#' There are `NUnrooted(n)` unrooted trees with _n_ leaves.
#' As such, each _n_-leaf tree can be uniquely identified by a non-negative
#' integer _x_ < `NUnrooted(n)`.
#'
#' This integer can be converted by a tree by treating it as a mixed-base
#' number, with bases 1, 3, 5, 7, &hellip; (2&nbsp;_n_ - 5).
#'
#' Each digit of this mixed base number corresponds to a leaf, and determines
#' the location on a growing tree to which that leaf should be added.
#'
#' We start with a two-leaf tree, and treat 0 as the origin of the tree.
#'
#' ```
#'
#' 0 ---- 1
#'
#' ```
#'
#' We add leaf 2 by breaking an edge and inserting a node (numbered
#' `2 + nTip - 1`).
#' In this example, we'll work up to a six-leaf tree; this node will be numbered
#' 2 + 6 - 1 = 7.
#' There is only one edge on which leaf 2 can be added.  Let's add node 7 and
#' leaf 2:
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
#' There are now three edges on which leaf 3 can be added.  Our options are:
#'
#' Option 0: the edge leading to 1;
#'
#' Option 1: the edge leading to 2;
#'
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
#' `1` is now the final digit of our mixed-base number.
#'
#' There are five places to add leaf 4:
#'
#' Option 0: the edge leading to 1;
#'
#' Option 1: the edge leading to 2;
#'
#' Option 2: the edge leading to 3;
#'
#' Option 3: the edge leading to 7;
#'
#' Option 4: the edge leading to 8.
#'
#' If we chose option 3, then `3` would be the penultimate digit of our
#' mixed-base number.
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
#' 2^64 unique 20-leaf trees.  As a `TreeNumber` is a 64-bit integer,
#' only trees with at most 19 leaves can be accommodated.
#'
#' @param x Integer identifying the tree (see details).
#' @param nTip Integer specifying number of leaves in the tree.
#' @template tipLabelsParam
#' @param \dots Additional parameters for consistency with S3 methods (unused).
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
#' @references
#' \insertAllCited{}
#' @encoding UTF-8
#' @seealso Describe the shape of a tree topology, independent of leaf labels:
#' [`TreeShape()`]
#' @family tree generation functions
#' @name TreeNumber
#

#' @rdname TreeNumber
#'
#' @return `as.TreeNumber()` returns an object of class `TreeNumber`,
#' which comprises a numeric vector, whose elements represent successive
#' nine-digit chunks of the decimal integer corresponding to the tree topology
#' (in big endian order).  The `TreeNumber` object has attributes
#' `nTip` and `tip.label`.
#' @export
as.TreeNumber <- function(x, ...) UseMethod('as.TreeNumber')

.TT_MAX_TIP <- 19L
# Calculate with:
# base <- cumprod(as.integer64(seq(1, by = 2, length.out = 16)))
.TT_BASE <- structure(c(
  4.94065645841247e-324, 1.48219693752374e-323, 7.4109846876187e-323, 
  5.18768928133309e-322, 4.66892035319978e-321, 5.13581238851976e-320, 
  6.67655610507569e-319, 1.00148341576135e-317, 1.7025218067943e-316, 
  3.23479143290917e-315, 6.79306200910926e-314, 1.56240426209513e-312, 
  3.90601065523782e-311, 1.05462287691421e-309, 3.05840634305121e-308, 
  7.87815130491404e-296), class = "integer64")

#' @rdname TreeNumber
#' @export
as.TreeNumber.phylo <- function(x, ...) {
  x <- RootTree(x, 1)
  edge <- x[["edge"]]
  nTip <- NTip(x)
  if (nTip > .TT_MAX_TIP) {
    warning("Trees with > 19 leaves not uniquely identified ",
            "by 64-bit TreeNumbers")
  }

  edge <- Postorder(edge)
  structure(.Int64(edge_to_num(edge[, 1], edge[, 2], nTip)),
            nTip = nTip,
            tip.label = TipLabels(x),
            class = c('TreeNumber', 'integer64'))
}

#' @rdname TreeNumber
#' @export
as.TreeNumber.multiPhylo <- function(x, ...) {
  lapply(x, as.TreeNumber.phylo)
}

#' @rdname TreeNumber
#' @export
as.TreeNumber.character <- function(x, nTip, tipLabels = TipLabels(nTip), ...) {
  structure(as.integer64(x),
            nTip = nTip,
            tip.label = tipLabels,
            class = c('TreeNumber', 'integer64'))
}

#' @rdname TreeNumber
#' @export
as.TreeNumber.TreeNumber <- function(x, ...) x

#' @export
#' @rdname TreeNumber
as.TreeNumber.MixedBase <- function(x, ...) {
  if (NTip(x) > .TT_MAX_TIP - 3L) {
    stop("Trees with > 19 leaves not uniquely identified by 64-bit TreeNumbers")
  }
  
  # Return:
  structure(sum(as.integer(x) * .TT_BASE[rev(seq_along(x))]),
            nTip = NTip(x),
            tip.label = TipLabels(x),
            class = c("TreeNumber", "integer64"))
}

#' @rdname TreeNumber
#' @export
as.MixedBase.TreeNumber <- function(x, ...) {
  nTip <- NTip(x)
  outLength <- nTip - 3L
  baseLength <- min(outLength, length(.TT_BASE))
  
  base <- .TT_BASE[rev(seq_len(baseLength))]
  ret <- integer(baseLength)
  n <- as.integer64(x)
  for (i in seq_len(baseLength)) {
    baseI <- base[i]
    ret[i] <- as.integer(n %/% baseI)
    n <- n %% baseI
  }
  
  # Return:
  structure(c(integer(outLength - baseLength), ret),
            nTip = NTip(x),
            tip.label = TipLabels(x),
            class = "MixedBase")
}

#' @rdname TreeNumber
#' @export
as.MixedBase.integer64 <- function(x, tipLabels = NULL, ...) {
  baseLength <- if (x > max(.TT_BASE)) {
    length(.TT_BASE)
  } else {
    which.max(.TT_BASE > x) - 1L
  }
  tipLabels <- TipLabels(tipLabels)
  if (is.null(tipLabels)) {
    tipLabels <- TipLabels(baseLength + 3L)
  }
  nTip <- length(tipLabels)
  outLength <- nTip - 3L
  if (outLength < baseLength) {
    stop("Number of tips too low; is tipLabels truncated?")
  }
  
  base <- .TT_BASE[rev(seq_len(baseLength))]
  ret <- integer(baseLength)
  
  for (i in seq_len(baseLength)) {
    baseI <- base[i]
    ret[i] <- as.integer(x %/% baseI)
    x <- x %% baseI
  }
  
  # Return:
  structure(c(integer(outLength - baseLength), ret),
            nTip = outLength + 3L,
            tip.label = tipLabels,
            binary = TRUE,
            class = "MixedBase")
}

#' @rdname TreeNumber
#' @export
as.MixedBase.numeric <- function(x, tipLabels = NULL, ...) {
  as.MixedBase.integer64(as.integer64(x), tipLabels = tipLabels, ...)
}
                                  
#' @rdname TreeNumber
#' @template MRS
#' @return `as.phylo.numeric()` returns a tree of class `phylo`.
#'
#' @examples
#' as.phylo(0:2, nTip = 6, tipLabels = letters[1:6])
#'
#' @importFrom ape as.phylo
#' @export
as.phylo.numeric <- function(x, nTip = attr(x, 'nTip'),
                              tipLabels = attr(x, 'tip.label'), ...) {
  if (is.null(nTip)) {
    if (is.null(tipLabels)) {
      stop("Either nTip or tipLabels must be specified.")
    } else {
      nTip <- length(tipLabels)
    }
  }
  if (is.null(tipLabels)) tipLabels <- paste0('t', seq_len(nTip))
  if (nTip == 1) {
    SingleTaxonTree(tipLabels)
  } else {
    if (length(x) > 1) {
      structure(lapply(x, as.phylo.numeric, nTip = nTip, tipLabels = tipLabels),
                tip.label = tipLabels, class = 'multiPhylo')
    } else {
      edge <- RenumberEdges(num_to_parent(.Int64.to.C(x), nTip),
                            seq_len(nTip + nTip - 2L))
      structure(list(edge = do.call(cbind, edge),
                     tip.label = tipLabels,
                     Nnode = nTip - 1L),
                order = 'preorder',
                class = 'phylo')
    }
  }
}

# Copied from as.phylo.numeric except if length > 1
#' @export
as.phylo.integer64 <- function(x, nTip = attr(x, 'nTip'),
                                tipLabels = attr(x, 'tip.label'), ...) {
  if (is.null(nTip)) {
    if (is.null(tipLabels)) {
      stop("Either nTip or tipLabels must be specified.")
    } else {
      nTip <- length(tipLabels)
    }
  }
  if (is.null(tipLabels)) tipLabels <- paste0('t', seq_len(nTip))
  if (nTip == 1) {
    SingleTaxonTree(tipLabels)
  } else {
    if (length(x) > 1) {
      ret <- vector('list', length(x))
      for (i in seq_along(x)) {
        ret[[i]] <- as.phylo(x[i], nTip = nTip, tipLabels = tipLabels)
      }
      ret <- structure(ret, tip.label = tipLabels, class = 'multiPhylo')
    } else {
      edge <- RenumberEdges(num_to_parent(.Int64.to.C(x), nTip),
                            seq_len(nTip + nTip - 2L))
      structure(list(edge = do.call(cbind, edge),
                     tip.label = tipLabels,
                     Nnode = nTip - 1L),
                order = 'preorder',
                class = 'phylo')
    }
  }
}

.Int64.to.C <- function(i64) {
  INT_MAX <- as.integer64(2147483647L)
  i64 <- as.integer64(i64)
  if (i64 > INT_MAX) {
    as.integer(c(i64 %/% INT_MAX, i64 %% INT_MAX))
  } else {
    as.integer(i64[1])
  }
}

#' @rdname TreeNumber
#' @export
as.phylo.TreeNumber <- function(x, nTip = attr(x, 'nTip'),
                                 tipLabels = attr(x, 'tip.label'), ...) {
  if (is.null(tipLabels)) tipLabels <- paste0('t', seq_len(nTip))
  edge <- RenumberEdges(num_to_parent(.Int64.to.C(x), nTip),
                        seq_len(nTip + nTip - 2L))
  structure(list(edge = do.call(cbind, edge),
                 tip.label = tipLabels,
                 Nnode = nTip - 1L),
            order = 'preorder',
            class = 'phylo')
}

#' @export
as.integer64.TreeNumber <- function(x, ...) {
  structure(x[1], class = 'integer64')
}

#' Print `TreeNumber` object
#'
#' S3 method for objects of class `TreeNumber`.
#'
#' @param x Object of class `TreeNumber`.
#' @param \dots Additional arguments for consistency with S3 method (unused).
#'
#' @export
print.TreeNumber <- function(x, ...) {
  nTip <- attr(x, 'nTip')
  cat("Phylogenetic tree number", .PrintedTreeNumber(x), "of",
      .PrintNUnrooted(nTip), "\n",
      nTip, "tips:", paste0(attr(x, 'tip.label')))
}

#' @keywords internal
.PrintNUnrooted <- function(n) {
  if (n < 15L || n > 19L) NUnrooted(n) else paste0(NUnrooted64(n))
}

# x: Object of class `TreeNumber`
#' @keywords internal
.PrintedTreeNumber <- function(x) {
  paste0(structure(x, class = 'integer64'))
}

#' @rdname TreeNumber
#' @export
as.MixedBase <- function(x, ...) UseMethod('as.MixedBase')

#' @rdname TreeNumber
#' @export
as.MixedBase.MixedBase <- function(x, ...) x

#' @rdname TreeNumber
#' @export
as.MixedBase.phylo <- function(x, ...) {
  x <- RootTree(x, 1)
  nTip <- NTip(x)
  edge <- x[["edge"]]
  
  edge <- Postorder(edge)
  structure(edge_to_mixed_base(edge[, 1], edge[, 2], nTip),
            nTip = nTip,
            tip.label = TipLabels(x),
            binary = TRUE,
            class = 'MixedBase')
}

#' @rdname TreeNumber
#' @export
as.MixedBase.multiPhylo <- function(x, ...) {
  lapply(x, as.MixedBase, ...)
}

#' @rdname TreeNumber
#' @export
as.phylo.MixedBase <- function(x, nTip = attr(x, 'nTip'),
                               tipLabels = attr(x, 'tip.label'), ...) {
  if (is.null(tipLabels)) {
    tipLabels <- paste0('t', seq_len(nTip))
  }
  edge <- RenumberEdges(mixed_base_to_parent(x, nTip),
                        seq_len(nTip + nTip - 2L))
  structure(list(edge = do.call(cbind, edge),
                 tip.label = tipLabels,
                 Nnode = nTip - 1L),
            order = 'preorder',
            class = 'phylo')
}