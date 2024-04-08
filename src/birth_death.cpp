#include <Rcpp/Light> // for min
#include "../inst/include/TreeTools/assert.h" /* for ASSERT */
#include "../inst/include/TreeTools/renumber_tree.h" /* for preorder */
#include <cmath> // for max, min
#include <random>
#include <vector>
#include <limits> // for infinity
#include <string> // for std::string&

using namespace Rcpp;

enum class Event {root, birth, /*death, survival*/ mutation, sampling};

struct bd_node {
  
  bd_node *child1;
  bd_node *child2;
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
    child1(nullptr), child2(nullptr), node_type(*node_type),
    event_type(event_type), time(time) {}
  
  void set_children(bd_node *child) {
    child1 = child;
    // child2 = nullptr;
  }
  
  void set_children(bd_node *first, bd_node *second) {
    child1 = first;
    child2 = second;
  }
  
  void sample() {
    event_type = Event::sampling;
  }
};

typedef bd_node* node_ptr;

inline void validate_dimension(const NumericMatrix &x, std::string x_name,
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

inline void validate_dimension(const NumericVector &x, std::string x_name,
                               const int *size) {
  if (x.size() != *size) {
    Rcpp::stop(x_name + " has length " + std::to_string(x.size()) +
      "; expecting " + std::to_string(*size));
  }
}

inline void validate_probability(const NumericVector &x, std::string x_name) {
  if (Rcpp::min(x) < 0) {
    Rcpp::stop(x_name + " contains entries < 0");
  }
  if (Rcpp::max(x) > 1) {
    Rcpp::stop(x_name + " contains entries > 1");
  }
}

inline void validate_sum_to_one(const NumericVector &x, std::string x_name) {
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

void process_node(node_ptr node, int* tip_i, int* node_i,
                  IntegerVector &edge1, IntegerVector& edge2,
                  NumericVector &length) {
  // tip_i is populated with the next UNUSED tip number.
  // Increment once used.
  // node_i is populated with the next UNUSED node number.
  // Decrement once used.
  
  switch (node->event_type) {
  
    case Event::root: {
      edge1.push_back((*node_i)--);
      node_ptr desc = node->child1;
      length.push_back(desc->time - node->time);
      process_node(desc, tip_i, node_i, edge1, edge2, length);
      break;
    }
    
    case Event::birth: {
      const int parent_node = (*node_i)--;
      const double parent_time = node->time;
      // 0. Terminate incoming edge
      edge2.push_back(parent_node);
      ASSERT(edge1.size() == edge2.size());
      ASSERT(edge1.size() == length.size());
      
      // 1. Create left edge and populate its parent
      edge1.push_back(parent_node);
      // 2. Set length of edge
      length.push_back(parent_time - node->child1->time);
      // 3. Process next event to continue left edge
      process_node(node->child1, tip_i, node_i, edge1, edge2, length);
      
      // 4. Do it all again for right edge
      edge1.push_back(parent_node);
      length.push_back(parent_time - node->child2->time);
      process_node(node->child2, tip_i, node_i, edge1, edge2, length);
      
      break;
    }
    case Event::mutation: {
      // 1. Add elapsed time to length of edge
      length[length.size() - 1] += (node->time - node->child1->time);
      // 2. Process next event to continue edge
      process_node(node->child1, tip_i, node_i, edge1, edge2, length);
      break;
    }
    case Event::sampling: {
      if (node->child1 == nullptr) {
        // Sampling marks end of lineage; label tip
        edge2.push_back((*tip_i)++);
        ASSERT(edge1.size() == edge2.size());
        ASSERT(edge1.size() == length.size());
      } else {
        // We have sampled an ancestor; create a zero-length edge and continue
        edge2.push_back((*node_i));
        ASSERT(edge1.size() == edge2.size());
        ASSERT(edge1.size() == length.size());
        
        edge1.push_back(*node_i);
        edge2.push_back((*tip_i)++);
        length.push_back(0.0);
        
        edge1.push_back((*node_i)--);
        length.push_back(node->time - node->child1->time);
        process_node(node->child1, tip_i, node_i, edge1, edge2, length);
      }
      break;
    }
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
    Rcpp::stop("`tMax` (" + std::to_string(tau) + ") must be non-negative");
  }
  
  
  std::mt19937 generator(rSeed[0]);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::exponential_distribution<double> exp1(1.0);
  
  // Reserve memory for "Set" 
  std::vector<std::vector<node_ptr>> set;
  set.reserve(n_types);
  
  for (int i = n_types; i--; ) {
    std::vector<node_ptr> inner_vector;
    inner_vector.reserve(nMax[0]);
    set.push_back(inner_vector);
  }
  
  // Draw RootNode with type a ~ pi
  double root_draw = uniform(generator);
  int root_type = n_types;
  while (root_draw >= 0) {
    root_draw -= pi[--root_type];
  }
  bd_node root_node(&root_type, Event::root, tau);
  set[root_type].push_back(&root_node);
  
  // Generate events until t < 0
  while (true) {
    
    // Get next event
    double next_offset = std::numeric_limits<double>::infinity();
    Event next_event = Event::root;
    int a = -1, b = -1;
    
    for (int i = n_types; i--; ) {
      const int n_of_type = set[i].size();
      if (!n_of_type) continue;
      // birth
      for (int j = n_types; j--; ) {
        const double offset = exp1(generator) / (n_of_type * lambda(i, j));
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
        const double offset = exp1(generator) / (n_of_type * gamma(i, j));
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
          // bd_node child(&i, Event::survival, 0); // Full model version
          
          // Under FE model, everything that survives is sampled
          bd_node child(&i, Event::sampling, 0);
          set[a][node]->set_children(&child);
          set[a][node] = &child;
        }
      }
      break;
    }
    
    // select parent node uniformly at random
    const int parent_i = set[a].size() * uniform(generator);
    switch (next_event) {
      case Event::birth: {
        bd_node 
          child1(&a, Event::birth, tau),
          child2(&b, Event::birth, tau)
        ;
        set[a][parent_i]->set_children(&child1, &child2);
        set[a][parent_i] = &child1;
        ASSERT(b >= 0);
        set[b].push_back(&child2);
        break;
      }
      /*case Event::death: {
        bd_node child(&a, Event::death, tau);
        set[a][parent_i].set_children(&child);
        set[a].erase(set[a].begin() + parent_i);
        break;
      }*/
      case Event::mutation: {
        bd_node child(&b, Event::mutation, tau);
        set[a][parent_i]->set_children(&child);
        set[a].erase(set[a].begin() + parent_i);
        ASSERT(b >= 0);
        set[b].push_back(&child);
        break;
      }
      case Event::sampling: {
        bd_node child(&a, Event::sampling, tau);
        set[a][parent_i]->set_children(&child);
        if (uniform(generator) < rA[a]) {
          // If sampled, it dies with probability ra
          set[a].erase(set[a].begin() + parent_i);
        } else {
          set[a][parent_i] = &child;
        }
        break;
      }
      // case Event::survival: {} // All survivors are sampled in FE model
      case Event::root: {
        Rcpp::stop("Event cannot be of type `root`. Please report this bug");
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
        node = set[type].erase(set[a].begin() + node);
      }
    }
  }*/
  // After the full phylogeny has been simulated, the tree is pruned to
  // remove all lineages and nodes that are not ancestral to a sampling event.
  // Moreover, birth events in which only one child lineage survives are
  // removed.
  // In the forward-equivalent simulation, pruning is not required.
  // prune_tree()
  
  // Now convert the nodes into a phylogenetic tree
  int tip = 1, node = -1;
  IntegerVector parent(0), child(0);
  NumericVector length(0);
  process_node(&root_node, &tip, &node, parent, child, length);
  
  ASSERT(parent.size() == child.size());
  ASSERT(parent.size() == length.size());
  
  for (int i = parent.size(); i--; ) {
    if (parent[i] < 0) {
      parent[i] = tip - parent[i];
    }
    if (child[i] < 0) {
      child[i] = tip - child[i];
    }
  }
  
  return TreeTools::preorder_weighted(parent, child, length);
}
