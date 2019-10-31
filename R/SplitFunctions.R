#' Subset of a split on fewer taxa
#'
#' @template splitsObjectParam
#' @param tips A vector specifying a subset of the tip labels applied to `split`.
#' @param keepAll logical specifying whether to keep entries that define trivial
#' splits (i.e. splits of zero or one tip) on the subset of tips.
#' @param unique logical specifying whether to remove duplicate splits.
#'
#' @return An object of class `Splits`, defined on `tips`.
#'
#' @examples
#'
#' splits <- as.Splits(PectinateTree(letters[1:9]))
#' efgh <- Subsplit(splits, tips = letters[5:8], keepAll = TRUE)
#' summary(efgh)
#'
#' TrivialSplits(efgh)
#'
#' Subsplit(splits, tips = letters[5:8], keepAll = FALSE)
#'
#'
#' @author Martin R. Smith
#'
#' @family split manipulation functions
#' @export
Subsplit <- function (splits, tips, keepAll = FALSE, unique = TRUE) {
  if (mode(splits) == 'list') {
    lapply(splits, Subsplit, tips = tips, keepAll = keepAll, unique = unique)
  } else if (length(splits) == 0) {
    ret <- splits
    attr(ret, 'nTip') <- length(tips)
    attr(ret, 'tip.label') <- if (mode(tips) == 'character') {
      tips
    } else {
      attr(splits, 'tip.label')[tips]
    }
    ret
  } else {
    allSplits <- as.Splits(as.logical(splits)[, tips])
    ret <- if (keepAll) allSplits else WithoutTrivialSplits(allSplits)

    # Return:
    if (unique) unique(ret) else ret
  }
}

#' Are splits trivial?
#'
#' @template splitsObjectParam
#'
#' @return Logical vector specifying whether each split in `splits` is trivial,
#' i.e. includes or excludes only a single tip or no tips at all.
#'
#' @author Martin R. Smith
#' @family split manipulation functions
#' @examples
#'
#' splits <- as.Splits(PectinateTree(letters[1:9]))
#' efgh <- Subsplit(splits, tips = letters[5:8], keepAll = TRUE)
#' summary(efgh)
#'
#' TrivialSplits(efgh)
#'
#' @export
TrivialSplits <- function (splits, nTip = attr(splits, 'nTip')) {
  inSplit <- TipsInSplits(splits, nTip)
  inSplit < 2L | inSplit > nTip - 2L
}

#' @describeIn TrivialSplits Remove trivial splits from a splits object
#' @export
WithoutTrivialSplits <- function (splits, nTip = attr(splits, 'nTip')) {
  splits[[!TrivialSplits(splits, nTip)]]
}

#' Which splits are compatible?
#'
#' @template splitsObjectParam
#' @param splits2 A second `Splits` object.
#'
#' @return A logical matrix specifying whether each split in `splits` is
#' compatible with each split in `splits2`.
#'
#' @examples
#'
#' splits <- as.Splits(BalancedTree(8))
#' splits2 <- as.Splits(PectinateTree(8))
#'
#' summary(splits)
#' summary(splits2)
#'
#' CompatibleSplits(splits, splits2)
#'
#' @author Martin R. Smith
#' @export
CompatibleSplits <- function (splits, splits2) {
  splits <- as.Splits(splits)
  nTip <- attr(splits, 'nTip')
  splits2 <- as.Splits(splits2, splits)
  apply(splits2, 1, function (split)
    apply(splits, 1, .CompatibleSplit, split, nTip))
}

#' @param a,b Integer representations of splits, from an row of a `Splits` object
#' @return logical vector stating whether splits are compatible
#' @describeIn CompatibleSplits Evaluate a single split pair
#' @keywords internal
#' @export
.CompatibleSplit <- function (a, b, nTip) {
  rawA <- .UintToRaw(a)
  rawB <- .UintToRaw(b)

  rawMask <- if (nTip %% 32) {
    .UintToRaw(c(2^(nTip %% 32) - 1, rep(2^32 - 1, nTip %/%32)))
  } else {
    .UintToRaw(rep(2^32 - 1, nTip %/%32))
  }
  .CompatibleRaws(rawA, rawB, rawMask)
}

#' @keywords internal
#' @export
.CompatibleRaws <- function (rawA, rawB, bitmask) {
  !any(as.logical(rawA & rawB)) ||
  !any(as.logical(rawA & !rawB)) ||
  !any(as.logical(!rawA & rawB)) ||
  !any(as.logical(!rawA & !rawB & bitmask))
}

.UintToRaw <- function (uint) {
  vapply(uint, function(x) as.raw(x %/% 2^c(24, 16, 8, 0) %% 256L), raw(4))
}

#' Probability of matching this well
#'
#' Calculates the probability that two random splits of the sizes provided
#' will be at least as similar as the two specified.
#'
#' @template split12Params
#'
#' @return The proportion of permissable informative splits
#' splitting the terminals into bipartitions of the sizes given,
#'  that match as well as `split1` and `split2` do.
#'
#' @examples
#' SplitMatchProbability(split1 = as.Splits(c(rep(TRUE, 4), rep(FALSE, 4))),
#'                       split2 = as.Splits(c(rep(TRUE, 3), rep(FALSE, 5))))
#'
#' @author Martin R. Smith
#' @export
SplitMatchProbability <- function (split1, split2) {

  if (Ntip(split1) != Ntip(split2)) stop("Splits pertain to different tips")
  if (!identical(attr(split1, 'tip.label'), attr(split2, 'tip.label'))) {
    stop("Sequence of tip labels must match")
  }

  split1 <- as.logical(split1)
  split2 <- as.logical(split2)
  partitions <- matrix(c(sum(split1 & split2),
                         sum(split1 & !split2),
                         sum(!split1 & split2),
                         sum(!split1 & !split2)
  ), 2, 2)
  #, dimnames=list(c('A1', 'B1'), c('A2', 'B2')))

  split1Size <- rowSums(partitions)
  split2Size <- colSums(partitions)
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
                         function (i) c(A1 - i, i, i + A2 - A1, B2 - i),
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
  sum(choices[ranking <= ranking[partitions[1, 2] + 1L - minA1B2]]) /
    choose(n, A1)
}

#' Distributions of tips consistent with a partition pair
#'
#' Number of terminal arrangements matching a specified configuration of
#' two partitions.
#'
#' Consider partitions that divide eight terminals, labelled A to H.
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
#'
#' @author Martin R. Smith
#' @export
NPartitionPairs <- function (configuration) {
  choose(sum(configuration[c(1, 3)]), configuration[1]) *
    choose(sum(configuration[c(2, 4)]), configuration[2])
}
#' @describeIn SplitMatchProbability The natural logarithm of the probability
LnSplitMatchProbability <- function(split1, split2) {
  log(SplitMatchProbability(split1, split2))
}
