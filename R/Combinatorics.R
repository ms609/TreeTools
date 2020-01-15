globalVariables(c('doubleFactorials', 'logDoubleFactorials'), 'TreeTools')

#' Double Factorial
#'
#' @param n Vector of integers.
#'
#' @return Returns the double factorial, n x (n - 2) x (n - 4) x (n - 6) x ...
#'
#' @examples {
#' DoubleFactorial (-4:0) # Return 1 if n < 2
#' DoubleFactorial (2) # 2
#' DoubleFactorial (5) # 1 x 3 x 5
#' exp(LnDoubleFactorial.int (8)) # 2 x 4 x 6 x 8
#'
#' }
#'
#' @template MRS
#' @family Double factorial
#' @export
DoubleFactorial <- function (n) {
  if (any(n > 300)) stop("301!! is too large to store as an integer. Use LnDoubleFactorial instead.")

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

# Memoizing this function makes it MUCH slower...
#' @describeIn DoubleFactorial Returns the logarithm of the double factorial.
#' @export
LnDoubleFactorial <- (function (n) {
  n[n < 2] <- 1 # Much faster than pmax
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

#' @rdname DoubleFactorial
#' @export
LogDoubleFactorial <- LnDoubleFactorial

#' @describeIn DoubleFactorial Slightly faster, when x is known to be length one
#' and below 50001
#' @export
LnDoubleFactorial.int <- function (n) {
  if (n < 2) {
    0
  } else {
    logDoubleFactorials[n]
  }
}

#' @rdname DoubleFactorial
#' @export
LogDoubleFactorial.int <- LnDoubleFactorial.int

#' Number of rooted/unrooted trees
#'
#' These functions return the number of rooted or unrooted trees consistent with
#' a given pattern of splits.
#'
#' Functions starting N return the number of rooted or unrooted trees, functions
#' starting Ln provide the natural logarithm of this number.
#' Calculations follow Carter _et al._ 1990, Theorem 2.
#'
#' @param tips Integer specifying the number of tips.
#' @param splits Integer vector listing the number of taxa in each tree
#' bipartition.
#'
#' @template MRS
#'
#' @references
#'  \insertRef{Carter1990}{TreeTools}
#'
#' @examples
#'   NRooted(10)
#'   NUnrooted(10)
#'   LnRooted(10)
#'   LnUnrooted(10)
#'   # Number of trees consistent with a character whose states are
#'   # 00000 11111 222
#'   NUnrootedMult(c(5,5,3))
#'
#' @export
NRooted     <- function (tips)  DoubleFactorial(tips + tips - 3L) # addition faster than 2*
#' @describeIn NRooted Number of unrooted trees
#' @export
NUnrooted   <- function (tips)  DoubleFactorial(tips + tips - 5L)
#' @describeIn NRooted  Log Number of unrooted trees
#' @export
LnUnrooted  <- function (tips) LnDoubleFactorial(tips + tips - 5L)
#' @describeIn NRooted  Log Number of unrooted trees (as integer)
#' @export
LnUnrooted.int <- function (tips) {
  ifelse(tips < 3L, 0, logDoubleFactorials[tips + tips - 5L])
}

#' @describeIn NRooted  Log Number of rooted trees
#' @export
LnRooted    <- function (tips) LnDoubleFactorial(tips + tips - 3L)
#' @describeIn NRooted  Log Number of rooted trees (as integer)
#' @export
LnRooted.int <- function (tips) {
  ifelse(tips < 2L, 0, logDoubleFactorials[tips + tips - 3L])
}

#' Number of trees one SPR step away
#'
#' `N1Spr` calculates the number of trees one subtree prune-and-regraft
#' operation away from a binary input tree using the formula given by Allen and
#' Steel (2001).
#'
#' `IC1Spr` calculates the information content of trees at this distance: i.e.
#' the entropy corresponding to the proportion of all possible _n_-tip trees
#' whose SPR distance is at most one from a specified tree.
#'
#' @param n Integer vector specifying the number of tips in a tree.
#'
#' @return `N1SPR` returns an integer vector.
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

#' @describeIn N1Spr Information content of trees 0 or 1 SPR step from tree
#'  with n tips.
#' @return `IC1SPR` returns an numeric vector.
#' @export
IC1Spr <- function(n) -log2((1L + N1Spr(n)) / NUnrooted(n))

#' @describeIn NRooted Log number of unrooted trees
#' @examples
#' LnUnrootedSplits(2,4)
#' LnUnrootedSplits(3,3)
#' @export
LnUnrootedSplits <- function (splits) {
  if ((nSplits <- length(splits)) < 2) return (LnUnrooted(splits));
  if (nSplits == 2) return (LnRooted(splits[1]) + LnRooted(splits[2]));
  return (LnUnrootedMult(splits))
}
#' @describeIn NRooted Number of unrooted trees consistent with a bipartition
#' split.
#' @examples
#' NUnrootedSplits(2,4)
#' NUnrootedSplits(3,3)
#' @family split information function
#' @export
NUnrootedSplits  <- function (splits) {
  if ((nSplits <- length(splits)) < 2) return (NUnrooted(splits));
  if (nSplits == 2) return (NRooted(splits[1]) * NRooted(splits[2]))
  return (NUnrootedMult(splits))
}
#' @describeIn NRooted Log unrooted mult
#' @export
LnUnrootedMult <- function (splits) {  # Carter et al. 1990, Theorem 2
  splits <- splits[splits > 0]
  totalTips <- sum(splits)

  # Return:
  LnDoubleFactorial(totalTips +  totalTips - 5L) -
    LnDoubleFactorial(2L * (totalTips - length(splits)) - 1L) +
    sum(LnDoubleFactorial(splits + splits - 3L))
}
#' @describeIn NRooted Number of unrooted trees consistent with a multi-partition
#' split
#' @export
NUnrootedMult  <- function (splits) {  # Carter et al. 1990, Theorem 2
  splits <- splits[splits > 0]
  totalTips <- sum(splits)

  # Return:
  round(DoubleFactorial(totalTips + totalTips - 5L) /
          DoubleFactorial(2L * (totalTips - length(splits)) - 1L)
        * prod(DoubleFactorial(splits + splits - 3L)))
}
