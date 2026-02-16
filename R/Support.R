#' Frequency of splits
#'
#' `SplitFrequency()` provides a simple way to count the number of times that
#' bipartition splits, as defined by a reference tree, occur in a forest of
#' trees. May be used to calculate edge ("node") support for majority consensus
#' or bootstrap trees.
#'
#' If multiple calculations are required, some time can be saved by using the
#' constituent functions (see examples).
#'
#' @param reference A tree of class `phylo`, a `Splits` object. If `NULL`,
#' the frequencies of all splits in `forest` will be returned.
#' @param forest A list of trees of class `phylo`, or a `multiPhylo` object; or a
#' `Splits` object. See
#' [vignette](https://ms609.github.io/TreeTools/articles/load-trees.html) for
#' possible methods of loading trees into R.
#'
#' @return `SplitFrequency()` returns the number of trees in `forest` that
#' contain each split in `reference`.
#' If `reference` is a tree of class `phylo`, then the sequence will correspond
#' to the order of nodes (use `ape::nodelabels()` to view).
#' Note that the three nodes at the root of the tree correspond to a single
#' split; see the example for how these might be plotted on a tree.
#'
#' @template exampleNodeSupport
#'
#' @template MRS
#' @family Splits operations
#' @export
SplitFrequency <- function(reference, forest) {
  if (is.null(reference)) {
    if (length(unique(lapply(lapply(forest, TipLabels), sort))) > 1) {
      stop("All trees must bear identical labels")
    }
    forestSplits <- do.call(c, as.Splits(forest, tipLabels = TipLabels(forest[[1]])))
    dup <- duplicated(forestSplits)
    ret <- forestSplits[[!dup]]
    logicals <- vapply(seq_along(forestSplits),
                       function(cf) ret %in% forestSplits[[cf]],
                       logical(sum(!dup)))
    count <- if (is.null(dim(logicals))) {
      sum(logicals)
    } else {
      rowSums(logicals)
    }
    attr(ret, "count") <- unname(count)
    ret
  } else {
    referenceSplits <- as.Splits(reference)
    refLabels <- attr(referenceSplits, "tip.label")
    forest <- lapply(forest, KeepTip, refLabels)
    forestSplits <- as.Splits(forest, tipLabels = refLabels)
  
    logicals <- vapply(forestSplits,
                       function(cf) referenceSplits %in% cf,
                       logical(length(referenceSplits)))
    ret <- if (is.null(dim(logicals))) {
      sum(logicals)
    } else {
      rowSums(logicals)
    }
    names(ret) <- rownames(referenceSplits)
  
    # Return:
    ret
  }
}

#' Label splits
#'
#' Labels the edges associated with each split on a plotted tree.
#'
#' As the two root edges of a rooted tree denote the same split, only the
#' rightmost (plotted at the bottom, by default) edge will be labelled.
#' If the position of the root is significant, add a tip at the root using
#' [`AddTip()`].
#'
#' @template treeParam
#' @param labels Named vector listing annotations for each split. Names
#' should correspond to the node associated with each split; see
#' [`as.Splits()`] for details.
#' If `NULL`, each splits will be labelled with its associated node.
#' @param unit Character specifying units of `labels`, if desired. Include a
#' leading space if necessary.
#' @param \dots Additional parameters to [`ape::edgelabels()`][ape::nodelabels].
#' @return `LabelSplits()` returns `invisible()`, after plotting `labels` on
#' each relevant edge of a plot (which should already have been produced using
#' `plot(tree)`).
#'
#' @examples
#' tree <- BalancedTree(LETTERS[1:5])
#' splits <- as.Splits(tree)
#' plot(tree)
#' LabelSplits(tree, as.character(splits), frame = "none", pos = 3L)
#' LabelSplits(tree, TipsInSplits(splits), unit = " tips", frame = "none",
#'             pos = 1L)
#'
#' @template exampleNodeSupport
#'
#' @seealso Calculate split support: [`SplitFrequency()`]
#'
#' Colour labels according to value: [`SupportColour()`]
#' @importFrom ape edgelabels
#' @importFrom stats setNames
#' @family Splits operations
#' @export
LabelSplits <- function(tree, labels = NULL, unit = "", ...) {
  splits <- as.Splits(tree)
  if (is.null(labels)) {
    splitNames <- names(splits)
    labels <- setNames(splitNames, splitNames)
  } else if (length(names(labels)) == 0) {
    if (length(splits) == length(labels)) {
      labels <- setNames(labels, names(splits))
    } else {
      stop("`labels` must bear the same names as `as.Splits(tree)`")
    }
  }
  
  if (length(setdiff(names(labels), names(splits)))) {
    warning("Label names do not correspond to splits in tree",
            immediate. = TRUE)
  }
  whichEdge <- match(as.integer(names(labels)), tree[["edge"]][, 2])
  edgelabels(paste0(labels, unit), edge = whichEdge, ...)
  # Return:
  invisible()
}

#' Colour for node support value
#'
#' Colour value with which to display node support.
#'
#' @param support A numeric vector of values in the range 0--1.
#' @param show1 Logical specifying whether to display values of 1.
#' A transparent white will be returned if `FALSE`.
#' @param scale 101-element vector listing colours in sequence. Defaults to
#' a diverging \acronym{HCL} scale.
#' @param outOfRange Colour to use if results are outside the range 0--1.
#' @return `SupportColour()` returns the appropriate value from `scale`,
#' or `outOfRange` if a value is outwith the valid range.
#' @examples
#' SupportColour((-1):4 / 4, show1 = FALSE)
#'
#' @template exampleNodeSupport
#'
#' @seealso Use in conjunction with [`LabelSplits()`] to colour split labels,
#' possibly calculated using [`SplitFrequency()`].
#'
#' @importFrom grDevices colorRampPalette
#' @export
SupportColour <- function(support,
                          show1 = TRUE,
                          # Equivalent to hcl.colors(101, "Blue-Red 2", rev = TRUE),
                          # but works safely on macOS
                          scale = colorRampPalette(c("#D33F6A", "#e2e2e2", "#4A6FE3"))(101),
                          outOfRange = "red") {
  sanitized <- support
  sanitized[!is.numeric(support) | support < 0 | support > 1] <- NA
  inRange <- !is.na(support) & support >= 0 & support <= 1 & support != ""
  ret <- character(length(support))
  ret[!inRange] <- outOfRange
  ret[inRange] <- scale[(sanitized[inRange] * 100) + 1L]
  if (!isTRUE(show1)) {
    ret[support == 1] <- "#ffffff00"
  }
  ret
}

#' @rdname SupportColour
#' @export
SupportColor <- SupportColour
