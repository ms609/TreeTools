#' Force a tree to match a constraint
#'
#' Modify a tree such that it matches a specified constraint.
#' This is at present a somewhat crude implementation that attempts to retain
#' much of the strucure of `tree` whilst guaranteeing compatibility with
#' each entry in `constraint`.
#'
#' @template treeParam
#' @template constraintParam
#'
#' @return `ImposeConstraint()` returns a tree of class `phylo`, consistent
#' with `constraint`.
#'
#' @examples
#' tips <- letters[1:9]
#' tree <- as.phylo(1, 9, tips)
#' plot(tree)
#'
#' constraint <- StringToPhyDat('0000?1111 000111111 0000??110', tips, FALSE)
#' plot(ImposeConstraint(tree, constraint))
#' @template MRS
#' @importFrom ape read.tree
#' @export
ImposeConstraint <- function (tree, constraint) {
  # This function is as efficient as it is elegant: i.e. not.
  # But it just about does the job.
  tree <- Preorder(tree)
  const <- AddUnconstrained(constraint,
                            setdiff(tree$tip.label, names(constraint)),
                            asPhyDat = FALSE)

  info <- apply(const, 2,
                function (x) SplitInformation(sum(x == '0'), sum(x == '1')))
  smallest <- ifelse(apply(const, 2, function (x) sum(x == '0') < sum(x == '1')),
                     '0', '1')

  tips <- tree$tip.label
  nTip <- length(tips)
  for (i in order(info)) {
    constI <- const[, i]
    collapsers <- constI == smallest[i]
    collapseNames <- names(collapsers[collapsers])
    if (length(collapseNames) < 2L) {
      warning("Could not apply constraint ", i,
              ". Is it incompatible or redundant?")
      next
    }
    collapsing <- apply(const[collapsers, , drop = FALSE], 2,
                        function (x) setdiff(x, '?')[1])

    const <- const[setdiff(rownames(const), collapseNames[-1]), , drop = FALSE]
    const[collapseNames[1], ] <- collapsing
    rownames(const)[match(collapseNames[1], rownames(const))] <- paste0(
      '(', paste0(collapseNames, collapse = ','), ')')

  }

  backbone <- Preorder(RenumberTips(read.tree(
    text = paste0('(', paste0(rownames(const), collapse = ','), ');')),
    tips))

  .ChildAtEnd <- function (x) {
    if (x <= nTip) x else .ChildAtEnd(edge[match(x, edge[, 1]), 2])
  }
  edge <- backbone$edge
  tomies <- table(edge[, 1], dnn = NULL)
  polytomies <- as.integer(names(tomies[tomies > 2]))
  for (node in polytomies) {
    nodeKids <- edge[edge[, 1] == node, 2]
    standIns <- vapply(nodeKids, .ChildAtEnd, 1)
    kept <- KeepTip(tree, standIns)$edge
    newNodes <- kept > length(standIns)
    kept[newNodes] <- kept[newNodes] - kept[1] + max(edge[, 1])
    kept[kept == kept[1]] <- node

    kept2 <- kept[, 2] # don't replace twice if standins[i - 1] < i
    for (i in seq_along(standIns)) {
      kept[, 2][kept2 == i] <- nodeKids[i]
    }
    edge <- rbind(edge[edge[, 1] != node, ], kept)
  }
  edge <- edge[order(edge[, 1]), ]
  backbone$edge <- RenumberTree(edge[, 1], edge[, 2])
  backbone$Nnode <- max(backbone$edge[, 1]) - nTip


  # Return:
  backbone
}

#' @describeIn ImposeConstraint Expand a constraint to include unconstrained
#' taxa.
#' @toAdd Character vector specifying taxa to add to constraint.
#' @asPhyDat Logical: if `TRUE`, return a `phyDat` object; if `FALSE`, return
#' a matrix.
#' @export
AddUnconstrained <- function (constraint, toAdd, asPhyDat = TRUE) {
  ret <- if (inherits(constraint, 'phyDat')) {
    PhyDatToMatrix(constraint)
  } else {
    constraint
  }
  ret <- rbind(ret, matrix('?', length(toAdd), dim(ret)[2],
                           dimnames = list(toAdd, NULL)))

  # Return:
  if (asPhyDat) {
    MatrixToPhyDat(ret)
  } else {
    ret
  }
}
