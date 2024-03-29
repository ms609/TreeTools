#' Simulate trees under the birth-death process
#' 
#' @param pi Vector of probabilities associated with each of _d_ types
#' @param lambda Square numeric matrix giving rate at which a lineage of type
#' _a_ gives birth (via cladogenesis) to an additional lineage of type _b_.
#' Entries on the diagonal give rate of (cladogenic) birth without mutation.
#' @param mu Numeric vector giving death rate for a lineage of each type
#' @param psi Numeric vector giving sampling rate for a lineage of each type
#' @param rA Numeric vector giving probability that a lineage of each type
#' dies after sampling
#' @param gamma Rate of anagenic mutation (i.e. change of type without a birth).
#' @param rAL Numeric vector giving probability of dying after a concerted
#' sampling event.
#' 
#' @param tMax time at which to start simulation (in units of time ago)
#' @param nMax Abort simulation when more than `nMax` lineages exist
#' 
#' Uses the forward-equivalent simulation described by
#' \insertRef{Celentano2024}{TreeTools}
#' https://arxiv.org/pdf/2402.17153v1.pdf
#' @references \insertAllCited{}
#' @template MRS
#' @export
BirthDeath <- function() {
  .FullBirthDeath(...)
}

.GetNextEvent <- function(
    popSize,
    tau,
    ...)

.FullBirthDeath <- function(
    pi,
    lambda,
    mu,
    psi,
    rA = double(length(pi)),
    gamma,
    rho,
    rAL = double(length(pi)),
    tMax,
    nMax
  ) {
  d <- length(pi)
  stopifnot(ncol(lambda) == d)
  stopifnot(nrow(lambda) == d)
  stopifnot(length(mu) == d)
  stopifnot(length(psi) == d)
  stopifnot(length(rA) == d)
  stopifnot(min(rA) >= 0)
  stopifnot(max(rA) <= 1)
  stopifnot(length(rAL) == d)
  stopifnot(min(rAL) >= 0)
  stopifnot(max(rAL) <= 1)
  stopifnot(length(gamma) == d)
  
  tau <- tMax
  popSize <- 
  eventTypes <- c("birth", "death", "mutation", "sampling")
}
