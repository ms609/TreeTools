#' Sort a list of phylogenetic trees
#'
#' Trees are sorted alphabetically by their Newick representation after
#' numbering tips alphabetically and renumbering tres in [`Preorder`].
#'
#' This method of sorting is inelegant, and will be replaced in a future
#' version by a method that is more principled
#' ([#84](https://github.com/ms609/TreeTools/issues/84)).
#'
#' @param x,decreasing,na.last,\dots As in [`sort()`].
#'
#' @examples
#' sort(as.phylo(0:5, 7))
#' @template MRS
#' @export
sort.multiPhylo <- function (x, decreasing = FALSE, na.last = NA, ...) {
  sortable <- vapply(x, .Sortable, character(1))
  x[order(sortable, na.last = na.last, decreasing = decreasing)]
}

#' @rdname sort.multiPhylo
#' @param e1,e2 `phylo` objects to be compared.
#' @export
`==.phylo` <- function (e1, e2) isTRUE(all.equal(e1, e2))

#' @rdname sort.multiPhylo
#' @export
`<.phylo` <- function (e1, e2) .Sortable(e1) < .Sortable(e2)

#' @rdname sort.multiPhylo
#' @export
`>.phylo` <- function (e1, e2) .Sortable(e1) > .Sortable(e2)


#' @importFrom ape write.tree
.Sortable <- function (tree) {
  write.tree(Preorder(RenumberTips(tree, sort(tree$tip.label))))
}

#
# xtfrm.multiPhylo <- function (x) {
#
# }
#
#
#
# ## sort.list
#
# function (x, partial = NULL, na.last = TRUE, decreasing = FALSE,
#           method = c("auto", "shell", "quick", "radix"))
# {
#   decreasing <- as.logical(decreasing)
#   if (is.null(partial) && is.numeric(x) && !is.object(x) &&
#       length(x) > 0) {
#     if (.Internal(sorted_fpass(x, decreasing, na.last)))
#       return(seq_along(x))
#   }
#   method <- match.arg(method)
#   if (method == "auto" && (is.numeric(x) || is.factor(x) ||
#                            is.logical(x) || (is.object(x) && !is.atomic(x))) &&
#       is.integer(length(x))) {
#     method <- "radix"
#   }
#   if (!is.null(partial))
#     .NotYetUsed("partial != NULL")
#   if (method == "quick") {
#     if (is.factor(x))
#       x <- as.integer(x)
#     if (is.numeric(x))
#       return(sort(x, na.last = na.last, decreasing = decreasing,
#                   method = "quick", index.return = TRUE)$ix)
#     else stop("method = \"quick\" is only for numeric 'x'")
#   }
#   if (is.na(na.last)) {
#     x <- x[!is.na(x)]
#     na.last <- TRUE
#   }
#   if (method == "radix") {
#     return(order(x, na.last = na.last, decreasing = decreasing,
#                  method = "radix"))
#   }
#   if (!is.atomic(x))
#     stop("'x' must be atomic for 'sort.list', method \"shell\" and \"quick\"\nHave you called 'sort' on a list?")
#   .Internal(order(na.last, decreasing, x))
# }
#
# ## rank
# ##
# function (x, na.last = TRUE, ties.method = c("average",
#                                              "first", "last", "random", "max",
#                                              "min"))
# {
#   nas <- is.na(x)
#   nm <- names(x)
#   ties.method <- match.arg(ties.method)
#   if (is.factor(x)) {
#     x <- as.integer(x)
#   }
#   x <- x[!nas]
#   y <- switch(ties.method, average = , min = ,
#               max = .Internal(rank(x, length(x), ties.method)),
#               first = sort.list(sort.list(x)),
#               last = sort.list(rev.default(sort.list(x, decreasing = TRUE))),
#               random = sort.list(order(x, stats::runif(sum(!nas)))))
#   if (!is.na(na.last) && any(nas)) {
#     yy <- NA
#     NAkeep <- (na.last == "keep")
#     if (NAkeep || na.last) {
#       yy[!nas] <- y
#       if (!NAkeep)
#         yy[nas] <- (length(y) + 1L):length(yy)
#     }
#     else {
#       len <- sum(nas)
#       yy[!nas] <- y + len
#       yy[nas] <- seq_len(len)
#     }
#     y <- yy
#     names(y) <- nm
#   }
#   else names(y) <- nm[!nas]
#   y
# }
# rank(x, ties.method = "min", na.last = "keep")
#
#
# ## xtfrm
# y <- if (is.numeric(x)) unclass(x) else {
#   as.vector(rank(x, ties.method = "min", na.last = "keep"))
# }
# if (!is.numeric(y) || ((length(y) != length(x)) && !inherits(x,
#                                                              "data.frame")))
#   stop("cannot xtfrm 'x'")
# y
# }
#
#
# ## order
# decreasing = TRUE
# z <- list(trees)
# decreasing <- as.logical(decreasing)
# if (length(z) == 1L && is.numeric(x <- z[[1L]]) && !is.object(x) &&
#     length(x) > 0) {
#   if (.Internal(sorted_fpass(x, decreasing, na.last)))
#     return(seq_along(x))
# }
# method <- match.arg(method)
# if (any(vapply(z, is.object, logical(1L)))) {
#   z <- lapply(z, function(x) if (is.object(x))
#     as.vector(xtfrm(x))
#     else x)
#   return(do.call("order", c(z, list(na.last = na.last,
#                                     decreasing = decreasing, method = method))))
# }
# if (method == "auto") {
#   useRadix <- all(vapply(z, function(x) {
#     (is.numeric(x) || is.factor(x) || is.logical(x)) &&
#       is.integer(length(x))
#   }, logical(1L)))
#   method <- if (useRadix)
#     "radix"
#   else "shell"
# }
# if (method != "radix" && !is.na(na.last)) {
#   return(.Internal(order(na.last, decreasing, ...)))
# }
# if (method == "radix") {
#   decreasing <- rep_len(as.logical(decreasing), length(z))
#   return(.Internal(radixsort(na.last, decreasing, FALSE,
#                              TRUE, ...)))
# }
# if (any(diff((l.z <- lengths(z)) != 0L)))
#   stop("argument lengths differ")
# na <- vapply(z, is.na, rep.int(NA, l.z[1L]))
# ok <- if (is.matrix(na))
#   rowSums(na) == 0L
# else !any(na)
# if (all(!ok))
#   return(integer())
# z[[1L]][!ok] <- NA
# ans <- do.call("order", c(z, list(decreasing = decreasing)))
# ans[ok[ans]]
# }
