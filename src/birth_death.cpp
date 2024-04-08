#include <Rcpp/Lightest>
#include <cmath> // for max, min
#include <random>
#include <vector>

using Rcpp;

enum class Event {root, birth, death, mutation, sampling, survival};

struct bd_node {
  
  using enum class Event;
  
  bd_node *child_1;
  bd_node *child_2;
  int node_type;
  Event event_type;
  double time;
  
  /*
  Equivalent to : 
  bd_node(int node_type, event event_type, double time) {
    this.node_type = node_type;
    this.event_type = event_type;
    this.time = time;
  }
  */
  
  bd_node(int node_type, Event event_type, double time) :
    node_type(node_type), event_type(event_type), time(time) {}
  
  void set_children(bd_node *child) {
    child_1 = child;
    // child_2 = nullptr;
  }
  
  void set_children(bd_node *first, bd_node *second) {
    child_1 = first;
    child_2 = second;
  }
};

void validate_dimension(const NumericMatrix &x, const char x_name,
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
void validate_dimension(const NumericVector &x, const char x_name,
                        const int *size) {
  if (x.size() != size) {
    Rcpp::stop(sprintf("%s has length %i; expecting %i",
                       x_name, x.length(), *size)
  }
}

void validate_probability(const NumericVector &x, const char x_name) {
  if (std::min(x) < 0) {
    Rcpp::stop(sprintf("%s contains entries < 0", x_name));
  }
  if (std::max(x) > 1) {
    Rcpp::stop(sprintf("%s contains entries > 1", x_name));
  }
})

void validate_sum_to_one(const NumericVector &x, const char x_name) {
  const int n = n.size();
  double sum = 0.0;
  for (int i = n; i--; ) {
    sum += x[i];
  }
  if (sum != 1) {
    Rcpp::stop(sprintf("Sum of %s should be one, not %d", x_name, sum);
  }
}
    
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
  validate_probability(&pi, "pi");
  validate_sum_to_one(&pi, "pi");
  validate_dimension(&lambda, "lambda", &n_types);
  validate_dimension(&mu, "mu", &n_types);
  validate_dimension(&psi, "psi", &n_types);
  validate_dimension(&rho, "rho", &n_types);
  validate_probability(&rho, "rho");
  validate_dimension(&rA, "rA", &n_types);
  validate_probability(&rA, "rA");
  validate_dimension(&rAL, "rAL", &n_types);
  validate_probability(&rAL, "rAL");
  validate_dimension(&gamma, "gamma", &n_types);
  
  double tau = tMax[0];
  
  std::mt19937 generator(std::random_device{}());
  std::uniform_real_distribution<double> uniform(0.0, 1.0);

  // Reserve memory for "Set" 
  std::vector<std::vector<bd_node>> set;
  set.reserve(n_types);
  
  for (int i = n_types; i--; ) {
    std::vector<bd_type> inner_vector;
    inner_vector.reserve(nMax[0]);
    set.push_back(inner_vector);
  }
  
  // Draw RootNode with type a ~ pi
  double root_draw = uniform(generator);
  for (int root_type = n_types; root_type--; ) {
    root_draw -= pi[root_type];
    if (root_draw < 0) {
      // Insert(RootNode, Sa)
      set[root_type].push_back(root_type, Event::root, tau);
      break;
    }
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
