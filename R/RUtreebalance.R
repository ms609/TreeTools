# This file is adapted from code by Rob Noble, doi:10.5281/zenodo.5873857
# https://github.com/robjohnnoble/RUtreebalance

# MIT License
# 
# Copyright (c) 2021 Rob Noble
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# Find the parent of a node.
# 
# Example:
# edges1 <- data.frame(Parent = c(1, 1, 1, 3, 3), Identity = 2:6)
# .MoveUp(edges1, 3)
.MoveUp <- function(edges, identity) {
  if (!(identity %in% edges[["Identity"]]) &&
      !(identity %in% edges[["Parent"]])) {
    stop("Invalid identity.")
  }
  parent <- edges[edges[["Identity"]] == identity, "Parent"]
  if (length(parent) == 0) {
    # if identity is the root then don't move
    identity
  } else if (is.factor(parent)) {
    levels(parent)[parent]
  } else {
    parent 
  }
}

# Find the root node of a tree.
# 
# Example:
# edges1 <- data.frame(Parent = c(1,1,1,3,3), Identity = 2:6)
# .FindRoot(edges1)
.FindRoot <- function(edges) {
  start <- edges[["Parent"]][[1]] # reasonable guess
  if (is.factor(start)) {
    start <- levels(start)[start]
  }
  repeat {
    if (.MoveUp(edges, start) == start) break
    start <- .MoveUp(edges, start)
  }
  return(start)
}

# Get a list of all subtree sizes via depth-first search.
# Optionally provide root node i and adjacency list Adj if known.
# 
# Example:
# tree1 <- data.frame(Parent = c(1, 1, 1, 1, 2, 3, 4), 
#                     Identity = 1:7, 
#                     Population = c(1, rep(5, 6)))
# .GetSubtreeSizes(tree1)
.GetSubtreeSizes <- function(
    tree,
    i = which(tree[["Identity"]] == .FindRoot(tree[, 1:2])),
    Adj = .GetAdjacency(tree),
    Col = NULL, Cumul = NULL, is_leaf =  NULL) {
  n <- length(tree[["Identity"]])
  has_pops <- "Population" %in% colnames(tree)
  
  uniqueID <- unique(tree[["Identity"]])
  if (is.null(Col)) {
    Col <- setNames(rep("w", n), uniqueID)
  }
  if (is.null(Cumul)) {
    Cumul <- setNames(rep(NA, n), uniqueID)
  }
  if (is.null(is_leaf)) {
    is_leaf <- setNames(logical(n), uniqueID)
  }
  if (is.null(Adj[[i]])) {
    is_leaf[[i]] <- TRUE
  }
  for (j in Adj[[i]]) {
    if (Col[[j]] == "w"){
      L <- .GetSubtreeSizes(tree, j, Adj, Col, Cumul, is_leaf)
      Col<- L[["colour"]]
      Cumul <- L[["cumulative"]]
      is_leaf <- L[["is_leaf"]]
    }
  }
  Col[[i]] <- "b"
  Cumul[[i]] <- if (has_pops) {
    tree[["Population"]][[i]] + sum(Cumul[Adj[[i]]])
  } else {
    ifelse(is_leaf[[i]], 1, 0) + sum(Cumul[Adj[[i]]])
  }
  # Return:
  list("colour" = Col,"cumulative" = Cumul, "is_leaf" = is_leaf)
}

# Get adjacency list of a tree.
# 
# Example:
# tree1 <- data.frame(Parent = c(1, 1, 1, 1, 2, 3, 4), 
#                     Identity = 1:7, 
#                     Population = c(1, rep(5, 6)))
# .GetAdjacency(tree1)
.GetAdjacency <- function(tree) {
  n <- length(tree[["Identity"]])
  Adj <- vector(mode = "list", length = n)
  for (i in seq_len(n)) {
    if(tree$Parent[[i]] != tree$Identity[[i]]) {
      p <- match(tree$Parent[[i]], tree$Identity)
      Adj[[p]] <- c(Adj[[p]], i)
    }
  }
  
  # Return:
  Adj
}

# Add a row to the edges list to represent the root node (if not already present)
.AddRootRow <- function(tree) {
  start <- setdiff(tree[["Parent"]], tree[["Identity"]])
  if (length(start) > 1) {
    stop("Input dataframe is missing one or more rows")
  }
  if (length(start) > 0) { # add row for root node
    root_row <- if ("Population" %in% colnames(tree)) {
      data.frame(Parent = start, Identity = start, Population = 0)
    } else {
      data.frame(Parent = start, Identity = start)
    }
    tree <- rbind(root_row, tree)
  }
  
  # Return:
  tree
}

#' Robust universal tree balance index
#' 
#' Calculate tree balance index \ifelse{html}{\out{J<sup>1</sup>}}{\eqn{J^1}}
#' (when `nonRootDominance = FALSE`) or
#' \ifelse{html}{\out{J<sup>1c</sup>}}{\eqn{J^{1c}}}
#' (when `nonRootDominance = TRUE`) from \insertCite{Lemant2022}{TreeTools}.
#' 
#' If population sizes are not provided, then the function assigns size 0 to
#' internal nodes, and size 1 to leaves.
#' 
#' @param tree Either an object of class 'phylo', or a dataframe with column 
#' names Parent, Identity and (optionally) Population.
#' The latter is similar to `tree$edge`, where tree is an object of class
#' 'phylo'; the differences are in class (data.frame versus matrix) and column
#' names.
#' The dataframe may (but isn't required to) include a row for the root node.
#' If population sizes are omitted then internal nodes will be assigned
#' population size zero and leaves will be assigned population size one.
#' @param q Numeric between zero and one specifying sensitivity to type
#' frequencies.  If `q < 1`, the \ifelse{html}{\out{J<sup>q</sup>}}{\eqn{J^q}}
#' index - based on generalized entropy - will be returned; see 
#' \insertCite{Lemant2022;textual}{TreeTools}, page 1223.
#' @param nonRootDominance Logical specifying whether to use non-root dominance
#' factor.
#' 
#' @returns `J1Index()` returns a numeric specifying the
#' \ifelse{html}{\out{J<sup>1</sup>}}{\eqn{J^1}} index of `tree`.
#' \eqn{J^1(T) = 1} for a perfectly balanced tree; 
#' \eqn{J^1(T) = 0} for a pectinate (linear / caterpillar) tree.
#' 
#' @examples
#' # Using phylo object as input:
#' phylo_tree <- read.tree(text="((a:0.1)A:0.5,(b1:0.2,b2:0.1)B:0.2);")
#' J1Index(phylo_tree)
#' phylo_tree2 <- read.tree(text='((A, B), ((C, D), (E, F)));')
#' J1Index(phylo_tree2)
#' 
#' # Using edges lists as input:
#' tree1 <- data.frame(Parent = c(1, 1, 1, 1, 2, 3, 4),
#'                     Identity = 1:7,
#'                     Population = c(1, rep(5, 6)))
#' J1Index(tree1)
#' tree2 <- data.frame(Parent = c(1, 1, 1, 1, 2, 3, 4),
#'                     Identity = 1:7,
#'                     Population = c(rep(0, 4), rep(1, 3)))
#' J1Index(tree2)
#' tree3 <- data.frame(Parent = c(1, 1, 1, 1, 2, 3, 4),
#'                     Identity = 1:7,
#'                     Population = c(0, rep(1, 3), rep(0, 3)))
#' J1Index(tree3)
#' cat_tree <- data.frame(Parent = c(1, 1:14, 1:15, 15),
#'                        Identity = 1:31,
#'                        Population = c(rep(0, 15), rep(1, 16)))
#' J1Index(cat_tree)
#'
#' # If population sizes are omitted then internal nodes are assigned population
#' # size zero and leaves are assigned population size one:
#' sym_tree1 <- data.frame(Parent = c(1, rep(1:15, each = 2)),
#'                        Identity = 1:31,
#'                        Population = c(rep(0, 15), rep(1, 16)))
#' # Equivalently:                        
#' sym_tree2 <- data.frame(Parent = c(1, rep(1:15, each = 2)),
#'                        Identity = 1:31)
#' J1Index(sym_tree1)
#' J1Index(sym_tree2)
#' @author Rob Noble, adapted by Martin R. Smith
#' @references \insertAllCited{}
#' @family tree characterization functions
#' @importFrom stats na.omit
#' @export
J1Index <- function(tree, q = 1, nonRootDominance = FALSE) {
  
  if (!is.na(tree)[[1]]) {
    if(inherits(tree, "phylo")) {
      # convert from phylo object to data frame
      tree <- as.data.frame(tree[["edge"]])
      colnames(tree) <- c("Parent", "Identity")
    }
    tree <- na.omit(tree) # remove any rows containing NA
    if(is.factor(tree[["Parent"]])) {
      tree[["Parent"]] <- levels(tree[["Parent"]])[tree[["Parent"]]]
    }
    if (is.factor(tree[["Identity"]])) {
      tree[["Identity"]] <- levels(tree[["Identity"]])[tree[["Identity"]]]
    }
    tree <- .AddRootRow(tree)
  }
  
  n <- length(tree[["Identity"]])
  if (n <= 1) {
    return(0)
  }
  
  # Adjacency list
  Adj <- .GetAdjacency(tree)
  
  # Get the list of all subtree sizes
  subtree_sizes <- .GetSubtreeSizes(tree, Adj = Adj)
  
  # Subtree sizes, including the root
  Cumul <- subtree_sizes[["cumulative"]]
  
  # Vector of internal nodes
  eff_int_nodes <- which(!subtree_sizes[["is_leaf"]])
  
  # Vector of leaves
  leaves <- which(subtree_sizes[["is_leaf"]])
  
  # If population sizes are missing then assign size 0 to internal nodes,
  # and size 1 to leaves:
  if(!("Population" %in% colnames(tree))) {
    tree[["Population"]] <- double(n)
    tree[["Population"]][leaves] <- 1
  }
  if (sum(tree$Population) <= 0) {
    stop("At least one node must have Population > 0")
  }
  J <- 0
  # Subtree sizes, excluding the root
  Star <- Cumul - tree[["Population"]]
  for (i in seq_len(n)){ # Loop over all nodes
    if (Star[[i]] > 0){ # Node has at least one child with non-zero size
      K <- 0
      if(length(Adj[[i]]) > 1) { # Otherwise i has only one child and its 
                                 # balance score is 0
        eff_children <- 0 # Number of children with non-zero size
        for (j in Adj[[i]]) {
          if (Cumul[[j]] > 0){ # Otherwise child j has a 0-sized subtree and 
                               # does not count
            eff_children <- eff_children + 1
            
            # p is the ratio of the child subtree size including the root
            # (root = the child) 
            # to the parent subtree size excluding the root
            p <- Cumul[[j]] / Star[[i]]
            
            # K is the sum of the node balance scores
            K <- K + if (q == 1) {-p * log(p)} else {p ^ q}
          }
        }
        # Non-root dominance factor:
        if (isTRUE(nonRootDominance)) {
          h_factor <- Star[[i]] / Cumul[[i]]
        } else {
          h_factor <- 1
        }
        
        # Normalize the sum of balance scores, adjust for non-root dominance, 
        # and then add the result to the index
        if (eff_children > 1) { # Exclude nodes that have only one child with 
                                # size greater than zero
          J_increment <- h_factor * Star[[i]] * if (q == 1) {
             K / log(eff_children)
          } else {
            (1 - K) * eff_children ^ (q - 1) / (eff_children ^ (q - 1) - 1)
          }
          J <- J + J_increment
        }
      }
    }
  }
  # normalize the index by dividing by the sum of all subtree sizes:
  if (length(eff_int_nodes) > 0) {
    J <- J / sum(Star[eff_int_nodes])
  }
  
  # Return:
  as.numeric(J)
}

#' @rdname J1Index
#' @export
JQIndex <- J1Index
