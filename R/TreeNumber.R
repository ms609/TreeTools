#' Unique integer indices for bifurcating tree topologies
#'
#' Functions converting between phylogenetic trees and their unique decimal
#' representation, based on a concept by John Tromp, employed in 
#' \insertCite{Li1996;textual}{TreeTools}.
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
#' `as.TreeNumber()` supports up to 51 leaves.
#' For trees with at most 19 leaves, the number fits in a 64-bit integer and
#' the `TreeNumber` is stored as an `integer64` (via the `bit64` package),
#' enabling arithmetic and exact round-tripping via `as.MixedBase()`.
#' For trees with 20–51 leaves, there are more than 2^64 distinct topologies,
#' so the tree number is stored as a decimal character string instead.
#'
#' Package developers can use the C++ header `TreeTools/tree_number.h`
#' (via `LinkingTo: TreeTools`) for the underlying 256-bit encoding
#' (`tree_num_t`) directly.
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
#' # Trees with 20–51 leaves are stored as decimal strings:
#' as.TreeNumber(BalancedTree(19))  # integer64-backed
#' as.TreeNumber(BalancedTree(51))  # character-backed
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
#' @family 'TreeNumber' utilities
#' @name TreeNumber
#

#' @rdname TreeNumber
#'
#' @return `as.TreeNumber()` returns an object of class `TreeNumber` with
#' attributes `nTip` and `tip.label`.
#' For trees with at most 19 leaves the underlying storage is a single
#' `integer64` value (class `c("TreeNumber", "integer64")`), enabling
#' `integer64` arithmetic and exact round-tripping through `as.MixedBase()`.
#' For trees with 20–51 leaves the number exceeds 2^64, so it is stored as a
#' decimal character string (class `c("TreeNumber", "character")`).
#' If `x` is a list of trees or a `multiPhylo` object,
#' `as.TreeNumber()` returns a corresponding list of `TreeNumber` objects.
#' @export
as.TreeNumber <- function(x, ...) UseMethod("as.TreeNumber")

# Maximum tips for integer64-backed TreeNumber (64-bit integers exhausted at 20 leaves).
.TT_MAX_TIP <- 19L
# Maximum tips supported by tree_number.h (256-bit integer, matches TREE_NUM_MAX_TIP).
.TT_MAX_TIP_FULL <- 51L

# Multiply a decimal string by a small positive integer, add a small
# non-negative integer, and return the result as a decimal string.
# Safe for multiplier and addend up to INT_MAX (~2.1e9): all intermediate
# values fit exactly in a 64-bit double (< 2^53).
.str_mul_add <- function(str, multiplier, addend) {
  m <- as.double(multiplier)
  digits <- as.integer(strsplit(str, "")[[1L]])
  n <- length(digits)
  # result_digits[k] holds the k-th least-significant output digit.
  result_digits <- integer(n + 12L)
  carry <- as.double(addend)
  for (k in seq_len(n)) {
    val <- digits[[n + 1L - k]] * m + carry
    result_digits[[k]] <- as.integer(val %% 10)
    carry <- floor(val / 10)
  }
  k <- n + 1L
  while (carry > 0) {
    result_digits[[k]] <- as.integer(carry %% 10)
    carry <- floor(carry / 10)
    k <- k + 1L
  }
  paste0(rev(result_digits[seq_len(k - 1L)]), collapse = "")
}

# Convert an INT_MAX-packed integer vector to a decimal string via Horner's
# method: result = chunks[1]; for each subsequent chunk c: result = result*INT_MAX + c.
.chunks_to_decimal <- function(chunks) {
  if (length(chunks) == 1L) return(as.character(chunks[[1L]]))
  result <- as.character(chunks[[1L]])
  for (i in seq_along(chunks)[-1L]) {
    result <- .str_mul_add(result, 2147483647L, chunks[[i]])
  }
  result
}

# Divide a decimal string by a small positive integer using long division.
# Returns list(quotient = character string, remainder = integer).
.str_divmod <- function(str, divisor) {
  divisor <- as.double(divisor)   # double avoids int32 overflow (max rem*10 ~ 2e10)
  digits <- as.integer(strsplit(str, "")[[1L]])
  quotient <- integer(length(digits))
  rem <- 0
  for (i in seq_along(digits)) {
    cur <- rem * 10 + digits[[i]]
    quotient[[i]] <- as.integer(cur %/% divisor)
    rem <- cur %% divisor
  }
  q_str <- sub("^0+", "", paste0(quotient, collapse = ""))
  if (nchar(q_str) == 0L) q_str <- "0"
  list(quotient = q_str, remainder = as.integer(rem))
}

# Convert a decimal string to an INT_MAX-packed integer vector.
# Inverse of .chunks_to_decimal(); used when converting character-backed
# TreeNumbers back to the form expected by num_to_parent().
.decimal_to_chunks <- function(str) {
  if (!nzchar(str) || str == "0") return(0L)
  chunks <- integer(0)
  while (str != "0") {
    dm <- .str_divmod(str, 2147483647L)
    chunks <- c(dm$remainder, chunks)
    str <- dm$quotient
  }
  if (length(chunks) == 0L) 0L else chunks
}

# Calculate with:
# base <- cumprod(as.integer64(seq(1, by = 2, length.out = 16)))
.TT_BASE <- as.integer64(c(
  "191898783962510625", "6190283353629375", "213458046676875", "7905853580625",
  "316234143225", "13749310575", "654729075", "34459425", "2027025", "135135",
  "10395", "945", "105", "15", "3", "1"))

.TTBases <- function (n) .TT_BASE[length(.TT_BASE) - n + seq_len(n)]

#' @export
as.TreeNumber.NULL <- function(x, ...) NULL

#' @rdname TreeNumber
#' @export
as.TreeNumber.phylo <- function(x, ...) {
  x <- RootTree(x, 1)
  edge <- x[["edge"]]
  nTip <- NTip(x)
  if (nTip > .TT_MAX_TIP_FULL) {
    stop("`as.TreeNumber()` supports at most ", .TT_MAX_TIP_FULL, " leaves.")
  }

  edge <- Postorder(edge)
  chunks <- edge_to_num(edge[, 1], edge[, 2], nTip)
  if (nTip <= .TT_MAX_TIP) {
    structure(.Int64(chunks),
              nTip = nTip,
              tip.label = TipLabels(x),
              class = c("TreeNumber", "integer64"))
  } else {
    structure(.chunks_to_decimal(chunks),
              nTip = nTip,
              tip.label = TipLabels(x),
              class = c("TreeNumber", "character"))
  }
}

#' @rdname TreeNumber
#' @export
as.TreeNumber.multiPhylo <- function(x, ...) {
  lapply(x, as.TreeNumber.phylo)
}

#' @rdname TreeNumber
#' @export
as.TreeNumber.character <- function(x, nTip, tipLabels = TipLabels(nTip), ...) {
  if (nTip <= .TT_MAX_TIP) {
    structure(as.integer64(x),
              nTip = nTip,
              tip.label = tipLabels,
              class = c("TreeNumber", "integer64"))
  } else {
    structure(x,
              nTip = nTip,
              tip.label = tipLabels,
              class = c("TreeNumber", "character"))
  }
}

#' @rdname TreeNumber
#' @export
as.TreeNumber.TreeNumber <- function(x, ...) x

#' @export
#' @rdname TreeNumber
as.TreeNumber.MixedBase <- function(x, ...) {
  if (NTip(x) > .TT_MAX_TIP) {
    stop("Converting a `MixedBase` to a `TreeNumber` requires `integer64` ",
         "arithmetic and is only supported for trees with <= ", .TT_MAX_TIP,
         " leaves. Use `as.phylo()` to convert larger `MixedBase` objects.")
  }

  # Return:
  structure(sum(as.integer(x) * .TTBases(length(x))),
            nTip = NTip(x),
            tip.label = TipLabels(x),
            class = c("TreeNumber", "integer64"))
}

#' @rdname TreeNumber
#' @export
as.MixedBase.TreeNumber <- function(x, ...) {
  nTip <- NTip(x)
  if (nTip > .TT_MAX_TIP) {
    stop("Converting a `TreeNumber` to `MixedBase` requires `integer64` ",
         "arithmetic and is only supported for trees with <= ", .TT_MAX_TIP,
         " leaves. Use `as.MixedBase(as.phylo(x))` for larger trees.")
  }
  outLength <- nTip - 3L
  baseLength <- min(outLength, length(.TT_BASE))

  base <- .TTBases(baseLength)
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
    which.max(rev(.TT_BASE) > x)
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
  
  base <- .TTBases(baseLength)
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
as.phylo.numeric <- function(x, nTip = attr(x, "nTip"),
                             tipLabels = attr(x, "tip.label"), ...) {
  if (is.null(nTip)) {
    if (is.null(tipLabels)) {
      stop("Either nTip or tipLabels must be specified.")
    } else {
      nTip <- length(tipLabels)
    }
  }
  if (nTip < 2) {
    if (nTip == 0) {
      ZeroTaxonTree()
    } else if (nTip == 1) {
      SingleTaxonTree(if (is.null(tipLabels)) TipLabels(nTip) else tipLabels)
    } else {
      stop("`nTip` may not be negative")
    }
  } else {
    if (is.null(tipLabels)) {
      tipLabels <- TipLabels(nTip)
    }
    if (length(x) > 1) {
      structure(lapply(x, as.phylo.numeric, nTip = nTip, tipLabels = tipLabels),
                tip.label = tipLabels, class = "multiPhylo")
    } else {
      edge <- RenumberEdges(num_to_parent(.Int64.to.C(x), nTip),
                            seq_len(nTip + nTip - 2L))
      .PreorderTree(
        edge = do.call(cbind, edge),
        tip.label = tipLabels,
        Nnode = nTip - 1L
      )
    }
  }
}

# Copied from as.phylo.numeric except if length > 1
#' @export
as.phylo.integer64 <- function(x, nTip = attr(x, "nTip"),
                                tipLabels = attr(x, "tip.label"), ...) {
  if (is.null(nTip)) {
    if (is.null(tipLabels)) {
      stop("Either nTip or tipLabels must be specified.")
    } else {
      nTip <- length(tipLabels)
    }
  }
  if (is.null(tipLabels)) {
    tipLabels <- paste0("t", seq_len(nTip))
  }
  if (nTip == 1) {
    SingleTaxonTree(tipLabels)
  } else {
    if (length(x) > 1) {
      ret <- vector("list", length(x))
      for (i in seq_along(x)) {
        ret[[i]] <- as.phylo(x[i], nTip = nTip, tipLabels = tipLabels)
      }
      ret <- structure(ret, tip.label = tipLabels, class = "multiPhylo")
    } else {
      edge <- RenumberEdges(num_to_parent(.Int64.to.C(x), nTip),
                            seq_len(nTip + nTip - 2L))
      .PreorderTree(edge = do.call(cbind, edge),
                    tip.label = tipLabels,
                    Nnode = nTip - 1L)
    }
  }
}

.Int64.to.C <- function(i64) {
  INT_MAX <- as.integer64(2147483647L)
  i64 <- as.integer64(i64)
  if (i64 > INT_MAX) {
    if (i64 > INT_MAX * INT_MAX) {
      stop("Number too large for 64-bit representation")
    }
    as.integer(c(i64 %/% INT_MAX, i64 %% INT_MAX))
  } else {
    as.integer(i64[1])
  }
}

#' @rdname TreeNumber
#' @export
as.phylo.TreeNumber <- function(x, nTip = attr(x, "nTip"),
                                 tipLabels = attr(x, "tip.label"), ...) {
  if (is.null(tipLabels)) tipLabels <- paste0("t", seq_len(nTip))
  chunks <- if (inherits(x, "integer64")) {
    .Int64.to.C(x)
  } else {
    .decimal_to_chunks(as.character(x))
  }
  edge <- RenumberEdges(num_to_parent(chunks, nTip),
                        seq_len(nTip + nTip - 2L))
  .PreorderTree(
    edge = do.call(cbind, edge),
    tip.label = tipLabels,
    Nnode = nTip - 1L
  )
}

#' @export
as.integer64.TreeNumber <- function(x, ...) {
  if (inherits(x, "integer64")) {
    structure(x[1], class = "integer64")
  } else {
    nTip <- attr(x, "nTip")
    if (nTip > .TT_MAX_TIP) {
      stop("Cannot convert a ", nTip, "-leaf TreeNumber to integer64: ",
           "trees with > ", .TT_MAX_TIP, " leaves require more than 64 bits.")
    }
    as.integer64(as.character(x))
  }
}

#' Is an object a `TreeNumber` object?
#' 
#' @param x R object.
#' @return `is.TreeNumber()` returns a logical vector of length one specifying
#' whether `x` inherits the class `"TreeNumber"`.
#' @template MRS
#' @examples
#' is.TreeNumber(FALSE) # FALSE 
#' is.TreeNumber(as.TreeNumber(BalancedTree(5))) # TRUE
#' @family 'TreeNumber' utilities
#' @export
is.TreeNumber <- function(x) inherits(x, "TreeNumber")

#' Print `TreeNumber` object
#'
#' S3 method for objects of class `TreeNumber`.
#'
#' @param x Object of class `TreeNumber`.
#' @param \dots Additional arguments for consistency with S3 method (unused).
#'
#' @family 'TreeNumber' utilities
#' @export
print.TreeNumber <- function(x, ...) {
  nTip <- attr(x, "nTip")
  cat("Phylogenetic tree number", .PrintedTreeNumber(x), "of",
      .PrintNUnrooted(nTip), "\n",
      nTip, "tips:", paste0(attr(x, "tip.label")))
}

#' @keywords internal
.PrintNUnrooted <- function(n) {
  if (n < 15L || n > 19L) NUnrooted(n) else paste0(NUnrooted64(n))
}

# x: Object of class `TreeNumber`
#' @keywords internal
.PrintedTreeNumber <- function(x) {
  if (inherits(x, "integer64")) {
    paste0(structure(x, class = "integer64"))
  } else {
    as.character(x)   # already a decimal string for character-backed TreeNumbers
  }
}

#' @rdname TreeNumber
#' @export
as.MixedBase <- function(x, ...) UseMethod("as.MixedBase")

#' @rdname TreeNumber
#' @export
as.MixedBase.MixedBase <- function(x, ...) x

#' @export
as.MixedBase.NULL <- function(x, ...) NULL

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
            class = "MixedBase")
}

#' @rdname TreeNumber
#' @export
as.MixedBase.multiPhylo <- function(x, ...) {
  lapply(x, as.MixedBase, ...)
}

#' @rdname TreeNumber
#' @export
as.phylo.MixedBase <- function(x, nTip = attr(x, "nTip"),
                               tipLabels = attr(x, "tip.label"), ...) {
  if (is.null(tipLabels)) {
    tipLabels <- paste0("t", seq_len(nTip))
  }
  edge <- RenumberEdges(mixed_base_to_parent(x, nTip),
                        seq_len(nTip + nTip - 2L))
  
  # Return:
  .PreorderTree(
    edge = do.call(cbind, edge),
    tip.label = tipLabels,
    Nnode = nTip - 1L
  )
}
