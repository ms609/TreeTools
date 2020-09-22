#' @importFrom utils globalVariables
utils::globalVariables(c('doubleFactorials',
                         'logDoubleFactorials',
                         'log2DoubleFactorials'),
                       'TreeTools')

#' Double factorial
#'
#' Calculate the double factorial of a number, or its logarithm.
#'
#' @param n Vector of integers.
#'
#' @return Returns the double factorial, _n_ * (_n_ - 2) * (_n_ - 4) *
#' (_n_ - 6) * ...
#'
#' @examples
#' DoubleFactorial (-4:0) # Return 1 if n < 2
#' DoubleFactorial (2) # 2
#' DoubleFactorial (5) # 1 * 3 * 5
#' exp(LnDoubleFactorial.int (8)) # log(2 * 4 * 6 * 8)
#' @template MRS
#' @family double factorials
#' @export
DoubleFactorial <- function (n) {
  if (any(n > 300)) stop("301!! is too large to represent. Use LnDoubleFactorial() instead.")

  n[n < 2] <- 1
  doubleFactorials[n]

  #
  #odds <- as.logical(x %% 2)
  #
  #oddX <- x[odds]
  #xPlusOneOverTwo <- (oddX + 1) / 2
  #evenX <- x[!odds]
  #xOverTwo <- evenX / 2
  #
  #ret <- integer(length(x))
  #ret[odds] <- gamma(oddX + 1L) / (gamma(xPlusOneOverTwo) * 2^(xPlusOneOverTwo - 1L))
  #ret[!odds] <- evenX * gamma(xOverTwo) * 2^(xOverTwo - 1L)
  #
  ## Return:
  #ret
}

#' @describeIn DoubleFactorial Returns the exact double factorial as a 64-bit
#' `integer64`, for `n` < 34.
#' @examples
#' DoubleFactorial64(31)
#' @export
DoubleFactorial64 <- function (n) {
  if (n < 2L) 1 else as.integer64(n * DoubleFactorial64(n - 2L))
}

# Memoizing this function makes it MUCH slower...
#' @describeIn DoubleFactorial Returns the logarithm of the double factorial.
#' @export
LnDoubleFactorial <- (function (n) {
  n[n < 2L] <- 1 # Much faster than pmax
  if (any(n > 49999L)) {

    odds <- as.logical(n %% 2)

    oddN <- n[odds]
    nPlusOneOverTwo <- (oddN + 1) / 2
    evenN <- n[!odds]
    nOverTwo <- evenN / 2

    ret <- integer(length(n))
    ret[odds] <- lgamma(oddN + 1L) - (lgamma(nPlusOneOverTwo) + (nPlusOneOverTwo - 1) * log(2))
    ret[!odds] <- log(evenN) + lgamma(nOverTwo) + (nOverTwo - 1) * log(2)

    # Return:
    ret

  } else {
    # Return from cache
    logDoubleFactorials[n]
  }
})

#' @describeIn DoubleFactorial Returns the logarithm of the double factorial.
#' @export
Log2DoubleFactorial <- (function (n) {
  n[n < 2L] <- 1 # Much faster than pmax
  if (any(n > 49999L)) {

    odds <- as.logical(n %% 2)

    oddN <- n[odds]
    nPlusOneOverTwo <- (oddN + 1) / 2
    evenN <- n[!odds]
    nOverTwo <- evenN / 2

    ret <- integer(length(n))
    ret[odds] <-  ((lgamma(oddN + 1L) - lgamma(nPlusOneOverTwo)) / log(2)) - nPlusOneOverTwo + 1
    ret[!odds] <- log2(evenN) + (lgamma(nOverTwo) / log(2)) + (nOverTwo - 1)

    # Return:
    ret

  } else {
    # Return from cache
    log2DoubleFactorials[n]
  }
})

#' @rdname DoubleFactorial
#' @export
LogDoubleFactorial <- LnDoubleFactorial

#' @describeIn DoubleFactorial Slightly faster, when x is known to be length one
#' and below 50001
#' @export
LnDoubleFactorial.int <- function (n) {
  if (n < 2L) {
    0
  } else {
    logDoubleFactorials[n]
  }
}

#' @rdname DoubleFactorial
#' @export
LogDoubleFactorial.int <- LnDoubleFactorial.int

#' Number of trees
#'
#' These functions return the number of rooted or unrooted binary trees
#' consistent with a given pattern of splits.
#'
#' Functions starting `N` return the number of rooted or unrooted trees.
#' Replace this initial `N` with `Ln` for the natural logarithm of this number;
#' or `Log2` for its base 2 logarithm.
#'
#' Calculations follow Cavalli-Sforza & Edwards (1967) and
#' Carter _et al._ 1990, Theorem 2.
#'
#' @param tips Integer specifying the number of leaves.
#' @param \dots Integer vector, or series of integers, listing the number of
#' leaves in each split.
#'
#' @template MRS
#'
#' @references
#'  \insertRef{Carter1990}{TreeTools}
#'
#'  \insertRef{CavalliSforza1967}{TreeTools}
#'
#' @examples
#' NRooted(10)
#' NUnrooted(10)
#' LnRooted(10)
#' LnUnrooted(10)
#' Log2Unrooted(10)
#' # Number of trees consistent with a character whose states are
#' # 00000 11111 222
#' NUnrootedMult(c(5,5,3))
#'
#' @family tree information functions
#' @export
NRooted     <- function (tips) DoubleFactorial(tips + tips - 3L) # addition faster than 2*

#' @describeIn NRooted Number of unrooted trees
#' @export
NUnrooted   <- function (tips) DoubleFactorial(tips + tips - 5L)

#' @describeIn NRooted Exact number of rooted trees as 64-bit integer
#' (13 < `nTip` < 19)
#' @export
NRooted64 <- function (tips) DoubleFactorial64(tips + tips - 3L)

#' @describeIn NRooted Exact number of unrooted trees as 64-bit integer
#' (14 < `nTip` < 20)
#' @examples
#' NUnrooted64(18)
#' @export
NUnrooted64 <- function (tips) DoubleFactorial64(tips + tips - 5L)

#' @describeIn NRooted  Log Number of unrooted trees
#' @export
LnUnrooted  <- function (tips) LnDoubleFactorial(tips + tips - 5L)

#' @describeIn NRooted  Log Number of unrooted trees (as integer)
#' @export
LnUnrooted.int <- function (tips) {
  ifelse(tips < 3L, 0, logDoubleFactorials[tips + tips - 5L])
}

#' @rdname NRooted
#' @export
Log2Unrooted  <- function (tips) Log2DoubleFactorial(tips + tips - 5L)

#' @rdname NRooted
#' @export
Log2Unrooted.int <- function (tips) {
  ifelse(tips < 3L, 0, log2DoubleFactorials[tips + tips - 5L])
}

#' @describeIn NRooted  Log Number of rooted trees
#' @export
LnRooted    <- function (tips) LnDoubleFactorial(tips + tips - 3L)
#' @describeIn NRooted  Log Number of rooted trees (as integer)
#' @export
LnRooted.int <- function (tips) {
  ifelse(tips < 2L, 0, logDoubleFactorials[tips + tips - 3L])
}
#' @rdname NRooted
#' @export
Log2Rooted    <- function (tips) Log2DoubleFactorial(tips + tips - 3L)
#' @rdname NRooted
#' @export
Log2Rooted.int <- function (tips) {
  ifelse(tips < 2L, 0, log2DoubleFactorials[tips + tips - 3L])
}

#' Number of trees one SPR step away
#'
#' `N1Spr()` calculates the number of trees one subtree prune-and-regraft
#' operation away from a binary input tree using the formula given by Allen and
#' Steel (2001); `IC1Spr()` calculates the information content of trees at this
#' distance: i.e. the entropy corresponding to the proportion of all possible
#' _n_-tip trees whose SPR distance is at most one from a specified tree.
#'
#' @param n Integer vector specifying the number of tips in a tree.
#'
#' @return `N1Spr()` returns an integer vector denoting the number of trees one
#' SPR rearrangement away from the input tree..
#'
#' @examples
#' N1Spr(4:6)
#' IC1Spr(5)
#'
#' @references
#'  \insertRef{Allen2001}{TreeTools}
#'
#' @export
N1Spr <- function (n) ifelse(n > 3L, (n + n - 6L) * (n + n - 7L), 0L)

#' @rdname N1Spr
#' @return `IC1Spr()` returns an numeric vector giving the phylogenetic
#' information content of trees 0 or 1 SPR rearrangement from an _n_-leaf tree,
#' in bits.
#' @export
IC1Spr <- function(n) Log2Unrooted(n) - log2(1L + N1Spr(n))

#' @rdname NRooted
#' @examples
#' LnUnrootedSplits(c(2,4))
#' LnUnrootedSplits(3, 3)
#' @export
LnUnrootedSplits <- function (...) {
  splits <- c(...)

  if ((nSplits <- length(splits)) < 2L) {
    LnUnrooted(splits)
  } else if (nSplits == 2L) {
    LnRooted(splits[1]) + LnRooted(splits[2])
  } else {
    LnUnrootedMult(splits)
  }
}

#' @rdname NRooted
#' @examples
#' Log2UnrootedSplits(c(2,4))
#' Log2UnrootedSplits(3, 3)
#' @export
Log2UnrootedSplits <- function (...) {
  splits <- c(...)

  if ((nSplits <- length(splits)) < 2L) {
    Log2Unrooted(splits)
  } else if (nSplits == 2L) {
    Log2Rooted(splits[1]) + Log2Rooted(splits[2])
  } else {
    Log2UnrootedMult(splits)
  }
}

#' @describeIn NRooted Number of unrooted trees consistent with a bipartition
#' split.
#' @examples
#' NUnrootedSplits(c(2,4))
#' NUnrootedSplits(3, 3)
#' @family split information function
#' @export
NUnrootedSplits  <- function (...) {
  splits <- c(...)
  if ((nSplits <- length(splits)) < 2L) {
    NUnrooted(splits)
  } else if (nSplits == 2L) {
    NRooted(splits[1]) * NRooted(splits[2])
  } else {
    NUnrootedMult(splits)
  }
}

#' @rdname NRooted
#' @export
LnUnrootedMult <- function (...) {  # Carter et al. 1990, Theorem 2
  splits <- c(...)
  splits <- splits[splits > 0]
  totalTips <- sum(splits)

  # Return:
  LnDoubleFactorial(totalTips +  totalTips - 5L) -
    LnDoubleFactorial(2L * (totalTips - length(splits)) - 1L) +
    sum(LnDoubleFactorial(splits + splits - 3L))
}

#' @rdname NRooted
#' @export
Log2UnrootedMult <- function (...) {  # Carter et al. 1990, Theorem 2
  splits <- c(...)
  splits <- splits[splits > 0]
  totalTips <- sum(splits)

  # Return:
  sum(Log2DoubleFactorial(totalTips +  totalTips - 5L),
      -Log2DoubleFactorial(2L * (totalTips - length(splits)) - 1L),
      Log2DoubleFactorial(splits + splits - 3L))
}

#' @describeIn NRooted Number of unrooted trees consistent with a multi-partition
#' split.
#' @export
NUnrootedMult  <- function (...) {  # Carter et al. 1990, Theorem 2
  splits <- c(...)
  splits <- splits[splits > 0]
  nSplits <- length(splits)
  if (nSplits < 1L) {
    # Return:
    0L
  } else if (nSplits == 1L) {
    # Return:
    NUnrooted(splits)
  } else if (nSplits == 2L) {
    prod(DoubleFactorial(splits + splits - 3L))
  } else {
    totalTips <- sum(splits)

    numerator <- totalTips + totalTips - 5L
    denominator <- 2L * (totalTips - length(splits)) + 1L

    # Return:
    prod(seq(numerator, denominator, -2L),
         DoubleFactorial(splits + splits - 3L))
  }
}
