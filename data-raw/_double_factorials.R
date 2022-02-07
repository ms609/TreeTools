# Memoizing this function makes it MUCH slower...
#' @describeIn DoubleFactorial Returns the logarithm of the double factorial.
LogDoubleFactorial <- (function(x) {
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
log2DoubleFactorials <- logDoubleFactorials / log(2)

DoubleFactorial <- function(x) {
  if (any(x > 300)) stop("301!! is too large to represent numerically. Use LogDoubleFactorial instead.")


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

doubleFactorials <- c(1L, 2L, 3L, 8L, 15L, 48L, 105L, 384L, 945L,
                      3840L, 10395L, 46080L, 135135L, 645120L, 2027025L,
                      10321920L, 34459425L, 185794560L, 654729075L,
                      # Use integers where possible to avoid rounding errors
                      exp(logDoubleFactorials[20:300]))
# Greater than 300 -> "Inf"

usethis::use_data(logDoubleFactorials, log2DoubleFactorials, doubleFactorials,
                  internal = TRUE, overwrite = TRUE, compress = 'xz')
