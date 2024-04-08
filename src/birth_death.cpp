#include <Rcpp/Lightest>
#include <cmath>
using Rcpp;

void validate_dimension(const NumericMatrix x, const char x_name,
                        const int *size) {
  if (x.ncol() != size) {
    Rcpp::stop(sprintf("%s has %i columns; expecting %i",
                       x_name, x.ncol(), *size)
  }
  if (x.nrow() != size) {
    Rcpp::stop(sprintf("%s has %i rows; expecting %s",
                       x_name, n.nrow(), *size)
  }
}
void validate_dimension(const NumericVector x, const char x_name,
                        const int *size) {
  if (x.length() != size) {
    Rcpp::stop(sprintf("%s has length %i; expecting %i",
                       x_name, x.length(), *size)
  }
}

void validate_probability(const NumericVector x, const char x_name) {
  if (std::min(x) < 0) {
    Rcpp::stop(sprintf("%s contains entries < 0", x_name));
  }
  if (std::max(x) > 1) {
    Rcpp::stop(sprintf("%s contains entries > 1", x_name));
  }
})

// [[Rcpp::export]]
List birth_death(
    const NumericVector pi,
    const NumericMatrix lambda,
    const NumericVector mu,
    const NumericVector psi,
    const NumericVector rA,
    const NumericMatrix gamma,
    const NumericVector rho,
    const NumericVector rAL,
    const NumericVector tMax,
    const NumericVector nMax
) {
  const int n_types = pi.length();
  // Check input is valid
  validate_dimension(lambda, "lambda", &n_types);
  validate_dimension(mu, "mu", &n_types);
  validate_dimension(psi, "psi", &n_types);
  validate_dimension(rho, "rho", &n_types);
  validate_probability(rho, "rho");
  validate_dimension(rA, "rA", &n_types);
  validate_probability(rA, "rA");
  validate_dimension(rAL, "rAL", &n_types);
  validate_probability(rAL, "rAL");
  validate_dimension(gamma, "gamma", &n_types);
  
  double tau = tMax[0];

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
