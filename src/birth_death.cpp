#include <Rcpp/Lightest>
#include "../inst/include/TreeTools/assert.h" /* for ASSERT */
#include <cmath> // for max, min
#include <random>
#include <vector>
#include <limits> // for infinity
#include <string> // for std::string&

using namespace Rcpp;

enum class Event {root, birth, /*death,*/ mutation, sampling, survival};

struct bd_node {
  
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
  
  bd_node(int* node_type, Event event_type, double time) :
    node_type(*node_type), event_type(event_type), time(time) {}
  
  void set_children(bd_node *child) {
    child_1 = child;
    // child_2 = nullptr;
  }
  
  void set_children(bd_node *first, bd_node *second) {
    child_1 = first;
    child_2 = second;
  }
  
  void sample() {
    event_type = Event::sampling;
  }
};

inline void validate_dimension(const NumericMatrix &x, std::string& x_name,
                               const int *size) {
  if (x.ncol() != *size) {
    Rcpp::stop(x_name + "has " + std::to_string(x.ncol()) +
      " columns; expecting " + std::to_string(*size));
  }
  if (x.nrow() != *size) {
    Rcpp::stop(x_name + "has " + std::to_string(x.nrow()) +
      " rows; expecting " + std::to_string(*size));
  }
}
inline void validate_dimension(const NumericVector &x, std::string& x_name,
                               const int *size) {
  if (x.size() != *size) {
    Rcpp::stop(x_name + " has length " + std::to_sting(x.size) +
      "; expecting " + std::to_string(*size));
  }
}

inline void validate_probability(const NumericVector &x, std::string& x_name) {
  if (min(x) < 0) {
    Rcpp::stop(x_name + " contains entries < 0");
  }
  if (max(x) > 1) {
    Rcpp::stop(x_name + " contains entries > 1");
  }
}

inline void validate_sum_to_one(const NumericVector &x, std::string& x_name) {
  const int n = x.size();
  double sum = 0.0;
  for (int i = n; i--; ) {
    sum += x[i];
  }
  if (sum != 1) {
    Rcpp::stop("Sum of " + x_name + " should be one, not " + 
      std::to_string(sum));
  }
}

// [[Rcpp::export]]
List birth_death(
    const NumericVector pi,
    const NumericMatrix lambda,
    const NumericVector psi,
    const NumericVector rA,
    const NumericMatrix gamma,
    const NumericVector rho,
    const NumericVector tMax,
    const IntegerVector nMax,
    const IntegerVector rSeed
) {
  const int n_types = pi.length();
  // Check input is valid
  validate_probability(pi, "pi");
  validate_sum_to_one(pi, "pi");
  validate_dimension(lambda, "lambda", &n_types);
  // validate_dimension(&mu, "mu", &n_types);
  validate_dimension(psi, "psi", &n_types);
  validate_dimension(rho, "rho", &n_types);
  validate_probability(rho, "rho");
  validate_dimension(rA, "rA", &n_types);
  validate_probability(rA, "rA");
  validate_dimension(gamma, "gamma", &n_types);
  
  const int n_max = nMax[0];
  if (n_max < 1) {
    Rcpp::stop("`nMax` (" + std::to_string(n_max) + " must be at least 1");
  }
  
  double tau = tMax[0];
  if (tau < 0) {
    Rcpp::stop("`tMax` (" + std::to_strung(tau) + ") must be non-negative");
  }
  
  
  std::mt19937 generator(rSeed[0]);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::exponential_distribution<double> exp1(1.0);

  
  inline int random_index(const std::vector<bd_node>& x) {
    static std::uniform_real_distribution<double> uniform_ref = uniform;
    static std::mt19937 generator_ref = generator;
    
    return x.size() * uniform_ref(generator_ref);
  }
  
  // Reserve memory for "Set" 
  std::vector<std::vector<&bd_node>> set;
  set.reserve(n_types);
  
  for (int i = n_types; i--; ) {
    std::vector<&bd_node> inner_vector;
    inner_vector.reserve(nMax[0]);
    set.push_back(inner_vector);
  }
  
  // Draw RootNode with type a ~ pi
  double root_draw = uniform(generator);
  int root_type;
  for (int root_type = n_types; root_type--; ) {
    root_draw -= pi[root_type];
    if (root_draw < 0) {
      // Insert(RootNode, Sa)
      break;
    }
  }
  const bd_node root_node(&root_type, Event::root, tau);
  set[root_type].push_back(&root_node);
  
  // Generate events until t < 0
  while (true) {
    
    // Get next event
    double next_offset = std::numeric_limits<double>::infinity();
    Event next_event;
    int a, b;
    
    for (int i = n_types; i--; ) {
      const int n_of_type = set[i].size();
      if (!n_of_type) continue;
      // birth
      for (int j = n_types; j--; ) {
        const double offset = exp1(generator) / (n_of_type * lambda[i, j]);
        if (offset < next_offset) {
          next_offset = offset;
          next_event = Event::birth;
          a = i;
          b = j;
        }
      }
      // death
      /* Not required for forward-equivalent
      const double death_offset = exp1(generator) / (n_of_type * mu[i]);
      if (death_offset < next_offset) {
        next_offset = death_offset;
        next_event = Event::death;
        a = i;
      } */
      // mutation
      for (int j = n_types; j--; ) {
        if (i == j) continue;
        const double offset = exp1(generator) / (n_of_type * gamma[i, j]);
        if (offset < next_offset) {
          next_offset = offset;
          next_event = Event::mutation;
          a = i;
          b = j;
        }
      }
      // sampling
      const double sample_offset = exp1(generator) / (n_of_type * psi[i]);
      if (sample_offset < next_offset) {
        next_offset = sample_offset;
        next_event = Event::sampling;
        a = i;
      }
    }
    
    // Increment time; if negative, mark survivors and exit while loop
    tau -= next_offset;
    
    if (tau < 0) {
      for (int i = n_types; i--; ) {
        for (int node = set[a].size(); node--; ) {
          bd_node child(&i, Event::survival, 0);
          set[a][node].set_children(&child);
          set[a][node] = &child;
        }
      }
      break;
    }
    const int parent_i = random_index(set[a]);
    switch (next_event) {
      case Event::birth: {
        bd_node 
          child1(&a, Event::birth, tau),
          child2(&b, Event::birth, tau)
        ;
        set[a][parent_i].set_children(&child1, &child2);
        set[a][parent_i] = &child1;
        set[b].push_back(&child2);
        break;
      }
      /*
      case Event::death:
        bd_node child(&a, Event::death, tau);
        set[a][parent_i].set_children(&child);
        set[a].erase(parent_i);
        break;*/
      case Event::mutation: {
        bd_node child(&b, Event::mutation, tau);
        set[a][parent_i].set_children(&child);
        set[a].erase(parent_i);
        set[b].push_back(&child);
        break;
      }
      case Event::sampling: {
        bd_node child(&a, Event::sampling, tau);
        set[a][parent_i].set_children(&child);
        if (uniform(generator) < rA[a]) {
          // If sampled, it dies with probability ra
          set[a].erase(parent_i);
        } else {
          set[a][parent_i] = &child;
        }
        break;
      }
    }
    int n_lineages = 0;
    for (int i = n_types; i--; ) {
      n_lineages += set[i].size();
    }
    ASSERT(n_lineages);
    if (!n_lineages) {
      // TODO return something appropriate
      Rcpp::stop("Empty");
    } else if (n_lineages > n_max) {
      // TODO return something appropriate
      Rcpp::stop("Capacity exceeded");
    }
  }
  // Now we have left the loop, any remaining nodes have survived to the present
  // and can be sampled (with p = rho[type])
  // Under the f-e process, rho = 1. 
  /*
  for (int type = n_types; type--; ) {
    for (auto node = set[type].begin(); node != set[type]end(); ) {
      if (uniform(generator) < ra[type]) {
        node->sample();
        ++node;
      } else {
        // Remove survivor from set (and update iterator)
        // Preempts pruning tree
        node = set[type].erase(node);
      }
    }
  }*/
  // After the full phylogeny has been simulated, the tree is pruned to
  // remove all lineages and nodes that are not ancestral to a sampling event.
  // Moreover, birth events in which only one child lineage survives are
  // removed.
  // In the forward-equivalent simulation, pruning is not required.
  // prune_tree()
  
  // Now can convert the nodes into a phylogenetic tree
  return List::create();
}
