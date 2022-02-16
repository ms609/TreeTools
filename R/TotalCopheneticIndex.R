#' Total Cophenetic Index
#'
#' `TotalCopheneticIndex()` calculates the total cophenetic index
#' \insertCite{Mir2013}{TreeTools} for any tree, a measure of its balance;
#' `TCIContext()` lists its possible values.
#'
#' The Total Cophenetic Index is a measure of tree balance -- i.e. whether
#' a (phylogenetic) tree comprises symmetric pairs of nodes, or has a pectinate
#' 'caterpillar' shape.
#' The index has a greater resolution power than Sackin's and Colless' indices,
#' and can be applied to trees that are not perfectly resolved.
#'
#' For a tree with _n_ leaves, the Total Cophenetic Index can take values of
#' 0 to `choose(n, 3)`.
#' The minimum value is higher for a perfectly resolved (i.e. dichotomous) tree
#' (see Lemma 14 of Mir _et al._ 2013).
#' Formulae to calculate the expected values under the Yule and Uniform models
#' of evolution are given in Theorems 17 and 23.
#'
#' Full details are provided by \insertCite{Mir2013;textual}{TreeTools}.
#'
#' @template xPhylo
#'
#' @return
#' `TotalCopheneticIndex()` returns an integer denoting the total cophenetic index.
#'
#' `TCIContext()` returns a data frame detailing the maximum and minimum value
#' obtainable for the Total Cophenetic Index for rooted binary trees with the
#' number of leaves of the given tree, and the expected value under the Yule
#' and Uniform models.
#' The variance of the expected value is given under the Yule model, but cannot
#' be obtained by calculation for the Uniform model.
#'
#' @seealso
#' `cophen.index()` in the package
#' '[CollessLike](https://github.com/LuciaRotger/CollessLike)'
#' provides an alternative implementation of this index and its predecessors.
#'
#' @references \insertAllCited{}
#'
#' @examples
#' # Balanced trees have the minimum index for a binary tree;
#' # Pectinate trees the maximum:
#' TCIContext(8)
#' TotalCopheneticIndex(PectinateTree(8))
#' TotalCopheneticIndex(BalancedTree(8))
#' TotalCopheneticIndex(StarTree(8))
#'
#'
#' # Examples from Mir et al. (2013):
#' tree12 <- ape::read.tree(text='(1, (2, (3, (4, 5))));')  #Fig. 4, tree 12
#' TotalCopheneticIndex(tree12) # 10
#' tree8  <- ape::read.tree(text='((1, 2, 3, 4), 5);')      #Fig. 4, tree 8
#' TotalCopheneticIndex(tree8)  # 6
#' TCIContext(tree8)
#' TCIContext(5L) # Context for a tree with 5 leaves.
#'
#' @family tree characterization functions
#'
#' @encoding UTF-8
#' @template MRS
#' @export
TotalCopheneticIndex <- function(x) UseMethod('TotalCopheneticIndex')

.Depth <- function(parent, child) {
  depth  <- integer(max(parent))
  for (i in seq_along(parent)) {
    depth[child[i]] <- depth[parent[i]] + 1L
  }

  # Return:
  depth
}

#' @importFrom fastmatch %fin%
#' @export
TotalCopheneticIndex.phylo <- function(x) {
  nTip   <- NTip(x)
  gc(T)
  edge   <- Preorder(x)[["edge"]]
  gc(T)
  parent <- edge[, 1]
  child  <- edge[, 2]
  depth  <- .Depth(parent, child)

  ancestors <- lapply(seq_len(nTip),
                      function(node) ListAncestors(parent, child, node))

  lca.depth <- vapply(seq_len(nTip), function(i) {
    vapply(seq_len(nTip), function(j) {
      anc.i <- ancestors[[i]]
      anc.j <- ancestors[[j]]
      lca <- max(anc.i[anc.i %fin% anc.j])

      # Return:
      depth[lca]
    }, integer(1L))
  }, integer(nTip))

  # Return:
  sum(lca.depth[upper.tri(lca.depth)])
}

#' @export
TotalCopheneticIndex.list <- function(x) {
  vapply(x, TotalCopheneticIndex, integer(1))
}

#' @export
TotalCopheneticIndex.multiPhylo <- TotalCopheneticIndex.list

#' @rdname TotalCopheneticIndex
#' @export
TCIContext <- function(x) UseMethod('TCIContext')

#' @export
TCIContext.phylo <- function(x) {
  TCIContext.numeric(NTip(x))
}

.MCI <- function(n) { # Lemma 14 in Mir er al 2013
  if (n < 3L) return(0L)
  halfN <- n / 2L
  topHalf <- ceiling(halfN)
  btmHalf <- floor(halfN)
  .MCI(topHalf) + .MCI(btmHalf) +
    choose(topHalf, 2L) + choose(btmHalf, 2L)
}

#' @rdname TotalCopheneticIndex
#' @export
TCIContext.numeric <- function(x) {
  H  <- function(n) sum(1 / (seq_len(n)))
  H2 <- function(n) sum(1 / (seq_len(n) ^ 2))

  maximum <- choose(x, 3L)
  minimum <- .MCI(x)

  # Theorem 17
  uniform.expected <- choose(x, 2) / 2L *
    ((DoubleFactorial((x + x) - 2L) / DoubleFactorial((x + x) - 3L)) - 2L)
  yule.expected    <- (x * (x + 1)) - (2 * x * H(x))
  yule.variance    <- ((1 / 12) * (x^4 - (10 * x^3) + (131 * x^2) - (2 * x))) -
    (4 * x^2 * H2(x)) - (6 * x * H(x))

  # Return:
  data.frame(maximum, minimum, uniform.expected, yule.expected, yule.variance)
}
