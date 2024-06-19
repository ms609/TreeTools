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
# edges1 <- data.frame(Parent = c(1,1,1,3,3), Identity = 2:6)
# .MoveUp(edges1, 3)
.MoveUp <- function(edges, identity) {
  if (!(identity %in% edges$Identity) & !(identity %in% edges$Parent)) {
    stop("Invalid identity.")
  }
  parent <- edges[edges$Identity == identity, "Parent"]
  if (length(parent) == 0) {
    # if identity is the root then don't move
    identity
  } else if(is.factor(parent)) {
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
  start <- edges$Parent[1] # reasonable guess
  if(is.factor(start)) start <- levels(start)[start]
  repeat {
    if(.MoveUp(edges, start) == start) break
    start <- .MoveUp(edges, start)
  }
  return(start)
}

# Get a list of all subtree sizes via depth-first search.
# Optionally provide root node i and adjacency list Adj if known.
# 
# Example:
# tree1 <- data.frame(Parent = c(1,1,1,1,2,3,4), 
#                     Identity = 1:7, 
#                     Population = c(1, rep(5, 6)))
# .GetSubtreeSizes(tree1)
.GetSubtreeSizes <- function(tree,i=NULL,Adj=NULL,Col=NULL,Cumul=NULL,is_leaf=NULL){
  n<-length(tree$Identity)
  has_pops <- FALSE
  if("Population" %in% colnames(tree)) has_pops <- TRUE
  if(is.null(Adj)) Adj <- .GetAdjacency(tree)
  if(is.null(i)) i <- which(tree$Identity == .FindRoot(tree[,1:2]))
  if(is.null(Col)) {
    Col <- rep("w",n)
    names(Col) <- unique(tree$Identity)
  }
  if(is.null(Cumul)) {
    Cumul <- rep(NA,n)
    names(Cumul) <- unique(tree$Identity)
  }
  if(is.null(is_leaf)) {
    is_leaf <- rep(FALSE, n)
    names(is_leaf) <- unique(tree$Identity)
  }
  if(is.null(Adj[[i]])) is_leaf[i] <- TRUE
  for (j in Adj[[i]]){
    if (Col[j] == "w"){
      L <- .GetSubtreeSizes(tree,j,Adj,Col,Cumul,is_leaf)
      Col<- L$colour
      Cumul <- L$cumulative
      is_leaf <- L$is_leaf
    }
  }
  Col[i] <- "b"
  if(has_pops) {
    Cumul[i] <- tree$Population[i] + sum(Cumul[Adj[[i]]])
  } else {
    Cumul[i] <- ifelse(is_leaf[i] == TRUE, 1, 0) + sum(Cumul[Adj[[i]]])
  }
  return(list("colour" = Col,"cumulative" = Cumul, "is_leaf" = is_leaf))
}

# Get adjacency list of a tree.
# 
# Example:
# tree1 <- data.frame(Parent = c(1,1,1,1,2,3,4), 
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
  start <- setdiff(tree$Parent, tree$Identity)
  if(length(start) > 1) stop("Input dataframe is missing one or more rows")
  if(length(start) > 0) { # add row for root node
    if("Population" %in% colnames(tree)) {
      root_row <- data.frame(Parent = start, Identity = start, Population = 0)
      message("Assigning Population = 0 to the root node")
    }
    else root_row <- data.frame(Parent = start, Identity = start)
    tree <- rbind(root_row, tree)
  }
  return(tree)
}

#' Robust universal tree balance index
#' 
#' Calculate tree balance index J^1 (when nonrootdomfactor = FALSE) or
#' J^{1c} (when nonrootdomfactor = TRUE) from \insertRef{Lemant2022}{TreeTools}.
#' If population sizes are missing then the function assigns
#' size 0 to internal nodes, and size 1 to leaves.
#' 
#' @param tree Either an object of class 'phylo', or a dataframe with column 
#' names Parent, Identity and (optionally) Population.
#' The latter is similar to `tree$edge`, where tree is a phylo object;
#' the differences are in class (data.frame versus matrix) and column names.
#' The dataframe may (but isn't required to) include a row for the root node.
#' If population sizes are omitted then internal nodes will be assigned
#' population size zero and leaves will be assigned population size one.
#' 
#' @examples
# Using phylo object as input:
# phylo_tree <- read.tree(text="((a:0.1)A:0.5,(b1:0.2,b2:0.1)B:0.2);")
# J1_index(phylo_tree)
# phylo_tree2 <- read.tree(text='((A, B), ((C, D), (E, F)));')
# J1_index(phylo_tree2)
# 
# Examples using edges lists as input:
# tree1 <- data.frame(Parent = c(1,1,1,1,2,3,4),
#                     Identity = 1:7,
#                     Population = c(1, rep(5, 6)))
# J1_index(tree1)
# tree2 <- data.frame(Parent = c(1,1,1,1,2,3,4),
#                     Identity = 1:7,
#                     Population = c(rep(0, 4), rep(1, 3)))
# J1_index(tree2)
# tree3 <- data.frame(Parent = c(1,1,1,1,2,3,4),
#                     Identity = 1:7,
#                     Population = c(0, rep(1, 3), rep(0, 3)))
# J1_index(tree3)
# cat_tree <- data.frame(Parent = c(1, 1:14, 1:15, 15),
#                        Identity = 1:31,
#                        Population = c(rep(0, 15), rep(1, 16)))
# J1_index(cat_tree)
#
# If population sizes are omitted then internal nodes are assigned population size zero
# and leaves are assigned population size one:
# sym_tree1 <- data.frame(Parent = c(1, rep(1:15, each = 2)),
#                        Identity = 1:31,
#                        Population = c(rep(0, 15), rep(1, 16)))
# sym_tree2 <- data.frame(Parent = c(1, rep(1:15, each = 2)),
#                        Identity = 1:31)
# all.equal(J1_index(sym_tree1), J1_index(sym_tree1))
#' @author Rob Noble, adapted by Martin R. Smith
#' @references \insertAllCited{}
#' @family tree characterization functions
#' @export
J1Index <- function(tree, q = 1, nonrootdomfactor = FALSE) {
  
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
  # adjacency list
  Adj <- .GetAdjacency(tree)
  # get the list of all subtree sizes
  subtree_sizes <- .GetSubtreeSizes(tree, Adj = Adj)
  # subtree sizes, including the root
  Cumul <- subtree_sizes[["cumulative"]]
  # vector of internal nodes
  eff_int_nodes <- which(!subtree_sizes[["is_leaf"]])
  # vector of leaves
  leaves <- which(subtree_sizes[["is_leaf"]])
  
  # if population sizes are missing then assign size 0 to internal nodes,
  # and size 1 to leaves:
  if(!("Population" %in% colnames(tree))) {
    tree[["Population"]] <- double(n)
    tree[["Population"]][leaves] <- 1
  }
  if (sum(tree$Population) <= 0) {
    stop("At least one node must have Population > 0")
  }
  J <- 0
  # subtree sizes, excluding the root
  Star <- Cumul - tree[["Population"]]
  for (i in 1:n){ # loop over all nodes
    if (Star[i] > 0){ # if node has at least one child with non-zero size
      K <- 0
      if(length(Adj[[i]])>1){ # otherwise i has only one child and its balance score is 0
        eff_children <- 0 # number of children with non-zero size
        for (j in Adj[[i]]){
          if (Cumul[j]>0){ # otherwise child j has a 0-sized subtree and does not count
            eff_children <- eff_children+1
            # p is the ratio of the child subtree size including the root (root = the child) 
            # to the parent subtree size excluding the root
            p <- Cumul[j]/Star[i]
            # K is the sum of the node balance scores
            if(q == 1) {
              K <- K + -p*log(p)
            } else {
              K <- K + p^q
            }
          }
        }
        # non-root dominance factor:
        if(nonrootdomfactor) {
          h_factor <- Star[i] / Cumul[i]
        } else {
          h_factor <- 1
        }
        # normalize the sum of balance scores, adjust for non-root dominance, 
        # and then add the result to the index
        if(eff_children > 1) { # exclude nodes that have only one child with size greater than zero
          if(q == 1) {
            J <- J + h_factor * Star[i] * K / log(eff_children)
          } else {
            J <- J + h_factor * Star[i] * (1 - K) * eff_children^(q - 1) / (eff_children^(q - 1) - 1)
          }
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

