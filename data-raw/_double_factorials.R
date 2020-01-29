# Memoizing this function makes it MUCH slower...
#' @describeIn DoubleFactorial Returns the logarithm of the double factorial.
LogDoubleFactorial <- (function (x) {
  x[x < 2L] <- 1L
  odds <- as.logical(x %% 2L)

  oddX <- x[odds]
  xPlusOneOverTwo <- (oddX + 1L) / 2L
  evenX <- x[!odds]
  xOverTwo <- evenX / 2L

  ret <- integer(length(x))
  ret[odds] <- lgamma(oddX + 1L) -
    (lgamma(xPlusOneOverTwo) + (xPlusOneOverTwo - 1L) * log(2L))
  ret[!odds] <- log(evenX) + lgamma(xOverTwo) + (xOverTwo - 1L) * log(2L)

  # Return:
  ret
})

logDoubleFactorials <- vapply(seq_len(50000), LogDoubleFactorial, double(1))

DoubleFactorial <- function (x) {
  if (any(x > 300)) stop("301!! is too large to store as an integer. Use LogDoubleFactorial instead.")


  x[x < 2L] <- 1L
  odds <- as.logical(x %% 2L)

  oddX <- x[odds]
  xPlusOneOverTwo <- (oddX + 1L) / 2L
  evenX <- x[!odds]
  xOverTwo <- evenX / 2L

  ret <- integer(length(x))
  ret[odds] <- gamma(oddX + 1L) / (gamma(xPlusOneOverTwo) * 2L^(xPlusOneOverTwo - 1L))
  ret[!odds] <- evenX * gamma(xOverTwo) * 2^(xOverTwo - 1L)

  # Return:
  ret
}

doubleFactorials <- exp(logDoubleFactorials[seq_len(300L)]) # Greater than 300 -> "Inf"
usethis::use_data(logDoubleFactorials, doubleFactorials, internal=TRUE, overwrite=TRUE)
