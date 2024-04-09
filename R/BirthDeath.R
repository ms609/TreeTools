#' Simulate trees under the birth-death process
#' 
#' @param pi Vector of probabilities associated with each of _d_ types
#' @param lambda Square numeric matrix giving rate at which a lineage of type
#' _a_ gives birth (via cladogenesis) to an additional lineage of type _b_.
#' Entries on the diagonal give rate of (cladogenic) birth without mutation.
#' @param mu Numeric vector giving death rate for a lineage of each type
#' @param psi Numeric vector giving sampling rate for a lineage of each type
#' @param r Numeric vector giving probability that a lineage of each type
#' dies after sampling. Defaults to certainty (`1`), such that ancestors
#' will not be sampled.
#' @param gamma Square numeric matrix giving rate of anagenic mutation 
#' (i.e. change of type without a birth)
#' from each type to each other. Diagonal is ignored (treated as zero).
#' @param rho Numeric vector giving probability that each type will be sampled
#' if it survives to the present.
# @param t Numeric vector specifying times of occurrence of concerted
# sampling events.  Not currently supported.
# @param rho Numeric matrix; each column gives sampling probability of each
# type at the concerted sampling event represented by each row.
# @param rL Numeric matrix giving probability of dying after being sampled
# during each concerted sampling event.
#' 
#' @param tMax time at which to start simulation (in units of time before present)
#' @param nMax Numeric; simulation will terminate when more than `nMax` lineages
#' exist, to avoid interminable simulation times.
#' @param seed Integer with which to seed mt19937 generator in C++.
#' Override the default to set a stated value for reproducible results.
#' @param times Vector giving time steps at which to solve differential
#' equations.  More steps give higher precision but require more
#' memory and take longer to initialize.
#' 
#' Uses the forward-equivalent simulation described by
#' \insertRef{Celentano2024}{TreeTools}
#' https://arxiv.org/pdf/2402.17153v1.pdf
#' @examples
#' lambda <- rbind(fit = c(fit = 1, unfit = 1), unfit = c(0.25, 0.25))
#' mu <- c(fit = 0.25, unfit = 0.25)
#' gamma <- cbind(fit = c(fit = 0, unfit = 0.5),
#'               unfit = c(fit = 0.5, unfit = 0))
#' 
#' BirthDeath(
#'  pi = c(fit = 0.5, unfit = 0.5),
#'  lambda = lambda,
#'  mu = mu,
#'  psi = c(fit = 0, unfit = 0),
#'  gamma = gamma,
#'  tMax = 500,
#'  nMax = 500
#'  )
#' 
#' @references \insertAllCited{}
#' @template MRS
#' @importFrom deSolve ode
#' @export
BirthDeath <- function(
    pi,
    lambda,
    mu,
    psi,
    gamma = matrix(0, length(pi), length(pi)),
    r = rep_len(1, length(pi)),
    rho = rep_len(1, length(pi)),
    tMax = 10,
    nMax = 1e5,
    seed = sample.int(.Machine[["integer.max"]], 1),
    times = seq(from = 0, to = tMax, length.out = 1001)
    ) {
  
  nTypes <- length(pi)
  yini  <- 1 - rho
  if (times[[1]] != 0) {
    stop("`times` must start at 0")
  }
  steps <- length(times)
  if (times[[steps]] != tMax) {
    stop("`times` must stop at `tMax` (", tMax, ")")
  }
  if (!all(times == cummax(times))) {
    stop("`times` must increase monotonically")
  }
  dEa_dt <- function(t, y, parms) {
    dy <-
      rowSums(vapply(seq_along(y),
                     function(b) lambda[, b] * ((y * y[b]) - y), y)) +
      mu * (1 - y) +
      psi * (0 - y) +
      rowSums(vapply(seq_along(y), function(b) gamma[, b] * (y[b] - y), y))
    list(dy)
  }
  bigE <- ode(y = yini, times = times, func = dEa_dt, parms = NULL)[, -1]
  oneMinusE <- 1 - bigE
  
  piBits <- pi * (1 - bigE[steps, ])
  piFE <- piBits / sum(piBits)
  
  # Careful with the dimensions:
  # Dimension 1 = Time step tau; dim 2 = a; dim 3 = b
  lambdaFE <- vapply(seq_len(nTypes), function(a) 
    vapply(seq_len(nTypes), function(b) 
      oneMinusE[, b] * lambda[a, b], double(steps)),
    oneMinusE)
  
  psiFE <- psi / oneMinusE
  
  rFE <- r + ((1 - r) * bigE)
  
  gammaFE <- vapply(seq_len(nTypes), function(a) 
    vapply(seq_len(nTypes), function(b) 
      (oneMinusE[, b] / oneMinusE[, a]) * (
        gamma[a, b] + 
          if(a == b) bigE[, a] * lambda[a, b] else 0
      ), double(steps)), oneMinusE)
  
  # rhoFE = 1
  
  birth_death(
    piFE,
    lambdaFE,
    psiFE,
    rFE,
    gammaFE,
    times,
    nMax,
    seed
  )
}



.FullBirthDeath <- function(
    pi,
    lambda,
    mu,
    psi,
    r = double(length(pi)),
    gamma,
    rho = rep_len(1, length(pi)),
    tMax,
    nMax
  ) {
  d <- length(pi)
  stopifnot(ncol(lambda) == d)
  stopifnot(nrow(lambda) == d)
  stopifnot(length(mu) == d)
  stopifnot(length(psi) == d)
  stopifnot(length(rho) == d)
  stopifnot(min(rho) >= 0)
  stopifnot(max(rho) <= 1)
  stopifnot(length(r) == d)
  stopifnot(min(r) >= 0)
  stopifnot(max(r) <= 1)
  stopifnot(length(rL) == d)
  stopifnot(min(rL) >= 0)
  stopifnot(max(rL) <= 1)
  stopifnot(ncol(gamma) == d)
  stopifnot(nrow(gamma) == d)
  
  tau <- tMax
  types <- names(pi)
  if (is.null(types)) {
    types <- seq_along(pi)
  }
  nTypes <- length(types)
  set <- setNames(vector("list", nTypes), types)
  set[[sample(types, 1, prob = pi)]] <- "birth" # TODO what does a new node look like?
  eventTypes <- c("birth", "death", "mutation", "sampling")
  
  .MatrixMin <- function(mat) {
    #index <- which.min(mat)
    #c(floor(index / nTypes), index %% nTypes + 1)
    which(mat == min(mat), arr.ind = TRUE, useNames = FALSE)
  }
  
  .NewNode <- function(type, event, time) {
    
  }
    
  while(tau >= 0) {
    
    # Get next event
    popSize <- lengths(set)
    births <- rexp(nTypes * nTypes) / (lambda * popSize)
    deaths <- rexp(nTypes) / (mu * popSize)
    mutations <- rexp(nTypes * nTypes) / (gamma * popSize)
    samplings <- rexp(nTypes) / (psi * popSize)
    nextTime <- vapply(list(births, deaths, mutations, samplings), min, double(1))
    nextEvent <- which.min(nextTime)
    ab <- switch(
      nextEvent, 
      .MatrixMin(births),
      c(which.min(deaths), NA_real_),
      .MatrixMin(mutations),
      c(which.min(samplings), NA_real_)
    )
    
    tau <- tau - nextTime[nextEvent]
    eventType <- eventTypes[nextEvent]
    eventA <- ab[[1]]
    eventB <- ab[[2]]
    
    if (tau < 0) {
      for (a in types) {
        for (node in set[[a]]) {
          child <- .NewNode(type = a, event = "survival", time = 0)
          .SetChildren(node, child)
          .Replace(set[[a]], node, child)
        }
      }
      break;
    }
    switch(
      eventType,
      "birth" = {
        child1 <- .NewNode(type = eventA, event = "birth", time = tau)
        child2 <- .NewNode(type = eventB, event = "birth", time = tau)
        parent <- .GetRandom(set[[eventA]])
        .SetChildren(parent, c(child1, child2))
        .Replace(set[[eventA]], parent, child1)
        .Insert(set[[eventB]], child2)
      }, "death" = {
        child <- .NewNode(type = eventA, event = "death", time = tau)
        parent <- .GetRandom(set[[eventA]])
        .SetChildren(parent, child)
        .Remove(set[[eventA]], parent)
      }, "mutation" = {
        child <- .NewNode(type = eventB, event = "mutation", time = tau)
        parent <- .GetRandom(set[[eventA]])
        .SetChildren(parent, child)
        .Remove(set[[eventA]], parent)
        .Remove(set[[eventB]], child)
      }, "sampling" = {
        child <- .NewNode(type = eventA, event = "sampling", time = tau)
        parent <- .GetRandom(set[[eventA]])
        .SetChildren(parent, child)
        .Replace(set[[eventA]], parent, child)
      }
    )
    n <- sum(lengths(set))
    if (n <= 0) {
      return(-Inf)
    }
    if (n > nMax) {
      return(Inf)
    }
    
    if (n <= nMax) { # Redundant check, if I've understood the return criteria
      for (a in types) {
        for (node in set[[a]]) {
          if (Bernoulli(rho[[a]]) == 1) {
            .SetEventType(node, event = "sampling")
          }
        }
      }
      .PruneTree()
    }
  }
}
