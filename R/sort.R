#' Sort a list of phylogenetic trees
#'
#' Trees are sorted by their [mixed base representation](MixedBase),
#' treating their leaves in the order of their labels (i.e. alphabetically,
#' if leaves are labelled with text).
#'
#' @param x,decreasing,na.last,\dots As in [`sort()`].
#'
#' @examples
#' sort(as.phylo(5:0, 7))
#' @template MRS
#' @export
sort.multiPhylo <- function(x, decreasing = FALSE, na.last = NA, ...) {
  maxTip <- max(NTip(x))
  sortable <- vapply(x, .Sortable, integer(maxTip - 3), maxTip)
  # Return:
  x[do.call(order, as.data.frame(t(sortable)), ...)]
}

#' @rdname sort.multiPhylo
#' @param e1,e2 Objects to be compared.
#' @export
`==.phylo` <- function(e1, e2) {
  isTRUE(all.equal(.Comparable(e1), .Comparable(e2)))
}

.Comparable <- function(tree) {
  as.MixedBase(RenumberTips(tree, sort(tree$tip.label)))
}

.Sortable <- function(tree, maxTip = NTip(tree)) {
  c(rep(-1L, maxTip - NTip(tree)),
    as.integer(.Comparable(tree)))
}

#' @rdname sort.multiPhylo
#' @export
`<.phylo` <- function(e1, e2) {
  .Comparable(e1) < .Comparable(e2)
}

#' @rdname sort.multiPhylo
#' @export
`>.phylo` <- function(e1, e2) `<.phylo`(e2, e1)


#' @rdname sort.multiPhylo
#' @export
`==.MixedBase` <- function(e1, e2) {
  isTRUE(all.equal(e1, e2))
}

#' @rdname sort.multiPhylo
#' @export
`<.MixedBase` <- function(e1, e2) {
  s1 <- as.integer(e1)
  s2 <- as.integer(e2)
  if (length(s1) == length(s2)) {
    s1Less <- s1 < s2
    if (any(s1Less)) {
      s1More <- s1 > s2
      if (any(s1More)) {
        which.max(s1Less) < which.max(s1More)
      } else {
        TRUE
      }
    } else {
      FALSE
    }
  } else {
    length(s1) < length(s2)
  }
}

#' @rdname sort.multiPhylo
#' @export
`>.MixedBase` <- function(e1, e2) e2 < e1
