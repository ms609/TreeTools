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
#' SplitMatchProbability(split1 = c(rep(TRUE, 4), rep(FALSE, 4)),
#'                       split2 = c(rep(TRUE, 3), rep(FALSE, 5)))
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
