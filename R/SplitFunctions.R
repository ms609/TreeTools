#' Subset of a split on fewer leaves
#'
#' `Subsplit()` removes leaves from a `Splits` object.
#'
#' @template splitsObjectParam
#' @param tips A vector specifying a subset of the leaf labels applied to `split`.
#' @param keepAll logical specifying whether to keep entries that define trivial
#' splits (i.e. splits of zero or one leaf) on the subset of leaves.
#' @param unique logical specifying whether to remove duplicate splits.
#'
#' @return `Subsplit()` returns an object of class `Splits`, defined on the
#' leaves `tips`.
#'
#' @examples
#' splits <- as.Splits(PectinateTree(letters[1:9]))
#' splits
#' efgh <- Subsplit(splits, tips = letters[5:8], keepAll = TRUE)
#' summary(efgh)
#'
#' TrivialSplits(efgh)
#'
#' summary(Subsplit(splits, tips = letters[5:8], keepAll = FALSE))
#' @template MRS
#' 
#' @seealso [`KeepTip()`] is a less flexible but faster equivalent.
#'
#' @family split manipulation functions
#' @export
Subsplit <- function(splits, tips, keepAll = FALSE, unique = TRUE) {
  if (is.list(splits)) {
    lapply(splits, Subsplit, tips = tips, keepAll = keepAll, unique = unique)
  } else if (length(splits) == 0) {
    ret <- splits
    attr(ret, "nTip") <- length(tips)
    attr(ret, "tip.label") <- if (is.character(tips)) {
      tips
    } else {
      attr(splits, "tip.label")[tips]
    }
    ret
  } else {
    allSplits <- as.Splits(as.logical(splits)[, tips])
    ret <- if (keepAll) allSplits else WithoutTrivialSplits(allSplits)

    # Return:
    if (unique) unique(ret) else ret
  }
}

#' Identify and remove trivial splits
#'
#' `TrivialSplits()` identifies trivial splits (which separate one or zero
#' leaves from all others); `WithoutTrivialSplits()` removes them from a
#' `Splits` object.
#'
#' @template splitsObjectParam
#' @template nTipParam
#'
#' @return `TrivialSplits()` returns a logical vector specifying whether each
#' split in `splits` is trivial, i.e. includes or excludes only a single tip or
#' no tips at all.
#'
#' @template MRS
#' @family split manipulation functions
#' @examples
#' splits <- as.Splits(PectinateTree(letters[1:9]))
#' efgh <- Subsplit(splits, tips = letters[5:8], keepAll = TRUE)
#' summary(efgh)
#'
#' TrivialSplits(efgh)
#' @export
TrivialSplits <- function(splits, nTip = attr(splits, "nTip")) {
  inSplit <- TipsInSplits(splits)
  inSplit < 2L | inSplit > nTip - 2L
}

#' @rdname TrivialSplits
#' @return `WithoutTrivialSplits()` returns a `Splits` object with trivial
#' splits removed.
#' @examples
#' summary(WithoutTrivialSplits(efgh))
#' @export
WithoutTrivialSplits <- function(splits, nTip = attr(splits, "nTip")) {
  splits[[!TrivialSplits(splits, nTip)]]
}

#' Which splits are compatible?
#'
#' @template splitsObjectParam
#' @param splits2 A second `Splits` object.
#'
#' @return `CompatibleSplits` returns a logical matrix specifying whether each
#' split in `splits` is compatible with each split in `splits2`.
#'
#' @examples
#' splits <- as.Splits(BalancedTree(8))
#' splits2 <- as.Splits(PectinateTree(8))
#'
#' summary(splits)
#' summary(splits2)
#'
#' CompatibleSplits(splits, splits2)
#'
#' @template MRS
#' @export
CompatibleSplits <- function(splits, splits2) {
  splits <- as.Splits(splits)
  nTip <- attr(splits, "nTip")
  splits2 <- as.Splits(splits2, splits)
  apply(splits2, 1, function(split)
    apply(splits, 1, .CompatibleSplit, split, nTip))
}

#' @param a,b [Raw][raw] representations of splits, from a row of a `Splits`
#' object.
#' @return `.CompatibleSplit` returns a logical vector stating whether splits
#' are compatible.
#' @rdname CompatibleSplits
#' @keywords internal
#' @export
.CompatibleSplit <- function(a, b, nTip) {
  rawMask <- if (nTip %% 8L) {
    as.raw(c(2 ^ (nTip %% 8L) - 1L, rep.int(255L, nTip %/% 8)))
  } else {
    rep.int(as.raw(255), nTip %/% 8L)
  }
  .CompatibleRaws(a, b, rawMask)
}

#' @rdname CompatibleSplits
#' @param rawA,rawB Raw representations of splits.
#' @param bitmask Raw masking bits that do not correspond to tips.
#'
#' @return `.CompatibleRaws` returns a logical vector specifying whether input
#' raws are compatible.
#' @keywords internal
#' @export
.CompatibleRaws <- function(rawA, rawB, bitmask) {
  !any(as.logical(rawA & rawB)) ||
  !any(as.logical(rawA & !rawB)) ||
  !any(as.logical(!rawA & rawB)) ||
  !any(as.logical(!rawA & !rawB & bitmask))
}

#' Probability of matching this well
#'
#' (`Ln`)`SplitMatchProbability()`calculates the probability that two random
#' splits of the sizes provided will be at least as similar as the two
#' specified.
#'
#' @template split12Params
#'
#' @return `SplitMatchProbability()` returns a numeric giving the proportion
#' of permissible non-trivial splits that divide the terminals into bipartitions
#' of the sizes given, that match as well as `split1` and `split2` do.
#'
#' @examples
#' split1 <- as.Splits(c(rep(TRUE, 4), rep(FALSE, 4)))
#' split2 <- as.Splits(c(rep(TRUE, 3), rep(FALSE, 5)))
#' SplitMatchProbability(split1, split2)
#' @family split information functions
#' @template MRS
#' @export
SplitMatchProbability <- function(split1, split2) {

  if (NTip(split1) != NTip(split2)) stop("Splits pertain to different tips")

  split1 <- as.logical(split1)
  split2 <- as.logical(as.Splits(split2, split1))
  partitions <- c(sum(split1 & split2), sum(split1 & !split2),
                  sum(!split1 & split2), sum(!split1 & !split2))
  #, dimnames=list(c("A1", "B1"), c("A2", "B2")))

  split1Size <- .rowSums(partitions, 2, 2)
  split2Size <- .colSums(partitions, 2, 2)
  A1 <- split1Size[1]
  A2 <- split2Size[1]
  B2 <- split2Size[2]
  n <- sum(partitions)

  minA1B2 <- max(0, A1 - A2)
  iMax <- min(A1, B2)
  minB1A2 <- B2 - iMax

  nArrangements <- iMax - minA1B2 + 1L
  # It turns out that this is called a confusion matrix, association matrix
  # or contingency table; see Meila 2007
  arrangements <- vapply(minA1B2:iMax,
                         function(i) c(A1 - i, i, i + A2 - A1, B2 - i),
                         double(4))

  #H <- function(p) -sum(p[p > 0] * log(p[p > 0]))
  #jointEntropies <- apply(arrangements / n, 2, H)
  #mutualInformation <- H(c(A1,n - A1) / n) + H(c(A2, B2) / n) - jointEntropies

  extraTipsInPerfectMatch.A1A2.B1B2 <- sum(arrangements[c(1, 4), nArrangements])
  extraTipsInPerfectMatch.A1B2.B1A2 <- sum(arrangements[c(2, 3), 1L])

  ranking <- if (extraTipsInPerfectMatch.A1A2.B1B2 >
                 extraTipsInPerfectMatch.A1B2.B1A2) {
    c(seq.int(1L, nArrangements, 2L),
      rev(seq.int(2L, nArrangements, 2L)))
  } else if (extraTipsInPerfectMatch.A1A2.B1B2 <
             extraTipsInPerfectMatch.A1B2.B1A2) {
    c(seq.int(2L, nArrangements, 2L),
      rev(seq.int(1L, nArrangements, 2L)))
  } else {
    c(seq.int(1L, nArrangements, 2L),
      rev(seq.int(2L, nArrangements, 2L)) - 1)
  }

  choices <- apply(arrangements, 2, NPartitionPairs)

  # Return:
  sum(choices[ranking <= ranking[partitions[3] + 1L - minA1B2]]) / choose(n, A1)
}

#' @rdname SplitMatchProbability
#' @return `LnSplitMatchProbability()` returns the natural logarithm of the
#' probability.
#' @examples
#' LnSplitMatchProbability(split1, split2)
#' @export
LnSplitMatchProbability <- function(split1, split2) {
  log(SplitMatchProbability(split1, split2))
}


#' Extract tip labels
#'
#' `TipLabels()` extracts labels from an object: for example, names of taxa in
#' a phylogenetic tree or data matrix.  `AllTipLabels()` extracts all labels,
#' where entries of a list of trees may pertain to different taxa.
#'
#' @param x An object of a supported class (see Usage section above).
#' @param single Logical specifying whether to report the labels for the first
#' object only (`TRUE`), or for each object in a list (`FALSE`).
#'
#' @return `TipLabels()` returns a character vector listing the tip labels
#' appropriate to `x`. If `x` is a single integer, this will be a vector
#' `t1`, `t2` ... `tx`, to match the default of \code{\link[ape]{rtree}()}.
#'
#' @examples
#' TipLabels(BalancedTree(letters[5:1]))
#' TipLabels(5)
#'
#' data("Lobo")
#' head(TipLabels(Lobo.phy))
#'
#' AllTipLabels(c(BalancedTree(4), PectinateTree(8)))
#'
#' @family tree properties
#' @template MRS
#' @export
TipLabels <- function(x, single = TRUE) UseMethod("TipLabels")

#' @rdname TipLabels
#' @export
TipLabels.default <- function(x, single = TRUE) {
  tla <- attr(x, "tip.label")
  # TODO when require R4.1, tla %||% names(x) %||% x
  if (length(tla) > 0) {
    tla
  } else {
    nom <- names(x)
    if (length(nom) == 0) {
      x
    } else {
      nom
    }
  }
}

#' @rdname TipLabels
#' @export
TipLabels.matrix <- function(x, single = TRUE) colnames(x)

#' @rdname TipLabels
#' @export
TipLabels.logical <- function(x, single = TRUE) TipLabels.numeric(length(x))

#' @rdname TipLabels
#' @export
TipLabels.phylo <- function(x, single = TRUE) x[["tip.label"]]

#' @rdname TipLabels
#' @export
TipLabels.phyDat <- function(x, single = TRUE) names(x)

#' @rdname TipLabels
#' @export
TipLabels.MixedBase <- TipLabels.default

#' @rdname TipLabels
#' @export
TipLabels.TreeNumber <- TipLabels.default

#' @rdname TipLabels
#' @family Splits operations
#' @export
TipLabels.Splits <- TipLabels.default

#' @rdname TipLabels
#' @export
TipLabels.list <- function(x, single = FALSE) {
  if (!is.null(attr(x, "tip.label"))) return(attr(x, "tip.label"))
  xTipLabel <- x[["tip.label"]]
  if (!is.null(xTipLabel)) {
    if (is.list(xTipLabel) && !is.null(xTipLabel[["tip.label"]])) {
      return(xTipLabel[["tip.label"]])
    } else {
      return(xTipLabel)
    }
  }
  .ListLabels(x, single, TipLabels)
}

#' @rdname TipLabels
#' @export
TipLabels.multiPhylo <- function(x, single = FALSE) {
  xTipLabel <- x[["tip.label"]]
  if (!is.null(xTipLabel)) {
    if (is.list(xTipLabel) && !is.null(xTipLabel[["tip.label"]])) {
      return(xTipLabel[["tip.label"]])
    } else {
      return(xTipLabel)
    }
  }
  if (single) {
    firstEntry <- x[[1]]
    if (!is.null(firstEntry[["tip.label"]])) {
      return(firstEntry[["tip.label"]])
    }
  } else {
    .ListLabels(x, single, TipLabels.phylo)
  }
}

#' @rdname TipLabels
#' @export
TipLabels.character <- function(x, single = TRUE) {
  if (is.null(attr(x, "tip.label"))) {
    NextMethod("TipLabels", as.character(x))
  } else {
    attr(x, "tip.label")
  }
}

#' @rdname TipLabels
#' @export
TipLabels.numeric <- function(x, single = TRUE) {
  if (length(x) == 1L) {
    if (x < 0) {
      stop("`x` may not be negative")
    }
    paste0(rep_len("t", x), seq_len(x))
  } else {
    NextMethod("TipLabels", as.character(x))
  }
}

#' @rdname TipLabels
#' @export
TipLabels.phyDat <- function(x, single = TRUE) names(x)

#' @rdname TipLabels
#' @export
AllTipLabels <- function(x) UseMethod("AllTipLabels")

#' @rdname TipLabels
#' @export
AllTipLabels.list <- function(x) {
  unique(unlist(lapply(x, TipLabels)))
}

#' @rdname TipLabels
#' @export
AllTipLabels.multiPhylo <- AllTipLabels.list

#' @rdname TipLabels
#' @export
AllTipLabels.phylo <- function(x) TipLabels.phylo(x)

#' @rdname TipLabels
#' @export
AllTipLabels.Splits <- function(x) TipLabels.Splits(x)

#' @rdname TipLabels
#' @export
AllTipLabels.TreeNumber <- function(x) TipLabels.TreeNumber(x)

#' @rdname TipLabels
#' @export
AllTipLabels.matrix <- function(x) TipLabels.matrix(x)

#' @keywords internal
.ListLabels <- function(x, single, Func) {
  if (length(x)) {
    if (single) {
      Func(x[[1]])
    } else {
      ret <- lapply(x, Func)
      uniqueRet <- unique(ret)
      if (length(uniqueRet) == 1) uniqueRet[[1]] else ret
    }
  } else {
    # else Return:
    NULL
  }
}

#' Distributions of tips consistent with a partition pair
#'
#' `NPartitionPairs()` calculates the number of terminal arrangements matching
#' a specified configuration of two splits.
#'
#' Consider splits that divide eight terminals, labelled A to H.
#'
#' \tabular{rcll}{
#'   Bipartition 1:\tab ABCD:EFGH\tab A1 = ABCD\tab B1 = EFGH \cr
#'   Bipartition 2:\tab ABE:CDFGH\tab A2 = ABE\tab B2 = CDFGH
#' }
#'
#' This can be represented by an association matrix:
#'
#' \tabular{rll}{
#'      \tab *A2* \tab *B2* \cr
#' *A1* \tab AB  \tab C   \cr
#' *B1* \tab E   \tab FGH
#' }
#'
#' The cells in this matrix contain 2, 1, 1 and 3 terminals respectively; this
#' four-element vector (`c(2, 1, 1, 3)`) is the `configuration` implied by
#' this pair of bipartition splits.
#'
#' @param configuration Integer vector of length four specifying the number of
#' terminals that occur in both
#'   (1) splits A1 and A2;
#'   (2) splits A1 and B2;
#'   (3) splits B1 and A2;
#'   (4) splits B1 and B2.
#'
#' @return The number of ways to distribute `sum(configuration)` taxa according
#'  to the specified pattern.
#'
#' @examples
#' NPartitionPairs(c(2, 1, 1, 3))
#' @template MRS
#' @export
NPartitionPairs <- function(configuration) {
  choose(sum(configuration[c(1, 3)]), configuration[1]) *
    choose(sum(configuration[c(2, 4)]), configuration[2])
}
