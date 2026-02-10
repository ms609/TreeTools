#' Split matching
#'
#' `match()` returns a vector of the positions of (first) matches of splits in
#' its first argument in its second.
#' `%in%` is a more intuitive interface as a binary operator, which returns
#' a logical vector indicating whether there is a match or not for each
#' split in its left operand.
#'
#' @param x,table Object of class `Splits`.
#' @param nomatch Integer value that will be used in place of `NA` in the case
#' where no match is found.
#' @param incomparables Ignored. (Included for consistency with generic.)
# @param incomparables A vector of values that cannot be matched. Any value in
# `x` matching a value in this vector is assigned the `nomatch` value.
# For historical reasons, `FALSE` is equivalent to `NULL`.
#'
#' @return `match()` returns an integer vector specifying the position in
#' `table` that matches each element in `x`, or `nomatch` if no match is found.
#'
#' @examples
#' splits1 <- as.Splits(BalancedTree(7))
#' splits2 <- as.Splits(PectinateTree(7))
#'
#' match(splits1, splits2)
# Turns base functions into S3 generics that can handle `Splits` objects
# (as well as `integer64`s).  This follows equivalent functions in the
# "\pkg{bit64}" package.
#
#' @seealso Corresponding base functions are documented in
#' [`match()`][base::match].
#'
#' @aliases match,Splits,Splits-method
#' @family Splits operations
#' @keywords methods
#' 
#' @rdname match.Splits
#' @export
setMethod("match",
          signature(x = "Splits", table = "Splits"),
          function(x, table, nomatch, incomparables) {
            if(missing("nomatch")) {
              nomatch <- NA_integer_
            }
            nomatch <- as.integer(nomatch)
            if (length(nomatch) != 1L) {
              nomatch <- NA_integer_
            }
            vapply(seq_along(x), function(i) {
              ret <- which(table %in% x[[i]])
              if (length(ret) == 0) ret <- nomatch
              ret
            }, integer(1))
          })

.in.Splits <- function(x, table) {
  duplicated(c(x, table), fromLast = TRUE)[seq_along(x)]
}

#' @rdname match.Splits
setGeneric("match")


#' @rdname match.Splits
#' @export
#' @keywords methods
setMethod("%in%",
          signature(x = "Splits", table = "Splits"),
          .in.Splits)


#' @importFrom methods setOldClass
#' @importFrom ape as.phylo
setOldClass(c("phylo", "multiPhylo"))

#' Tree matching
#'
#' `match()` returns a vector of the positions of (first) matches of trees in
#' its first argument in its second.
#' `%in%` is a more intuitive interface as a binary operator, which returns
#' a logical vector indicating whether there is a match or not for each
#' tree in its left operand.
#'
#' @param x,table Object of class `phylo` or `multiPhylo`.
#' @param nomatch Integer value that will be used in place of `NA` in the case
#' where no match is found.
#' @param incomparables Ignored. (Included for consistency with generic.)
# @param incomparables A vector of values that cannot be matched. Any value in
# `x` matching a value in this vector is assigned the `nomatch` value.
# For historical reasons, `FALSE` is equivalent to `NULL`.
#'
#' @return `match()` returns an integer vector specifying the position in
#' `table` that matches each element in `x`, or `nomatch` if no match is found.
#'
#' @examples
#' tree1 <- BalancedTree(7)
#' trees <- c(PectinateTree(7), BalancedTree(7))
#'
#' match(tree1, trees)
#
#' @seealso Corresponding base functions are documented in
#' [`match()`][base::match].
#'
#' @family utility functions
#' @keywords methods
#' 
#' @rdname match.multiPhylo
#' @aliases match,phylo,phylo-method
#' @export
setMethod("match",
          signature(x = "phylo", table = "phylo"),
          function(x, table, nomatch = NA_integer_, incomparables = NULL) {
            if (all.equal(x, table)) {
              1L
            } else {
              if (missing("nomatch")) {
                NA_integer_
              } else {
                nomatch <- as.integer(nomatch)
                if (length(nomatch) != 1L) {
                  NA_integer_
                } else {
                  nomatch
                }
              }
            }
          })


#' @rdname match.multiPhylo
#' @export
setMethod("match",
          signature(x = "multiPhylo", table = "phylo"),
          function(x, table, nomatch = NA_integer_, incomparables = NULL) {
            if(missing("nomatch")) {
              nomatch <- NA_integer_
            }
            nomatch <- as.integer(nomatch[[1]])
            if (length(nomatch) != 1L) {
              nomatch <- NA_integer_
            }
            ifelse(x %in% table, 1L, nomatch)
          })

#' @rdname match.multiPhylo
#' @export
setMethod("match",
          signature(x = "phylo", table = "multiPhylo"),
          function(x, table, nomatch = NA_integer_, incomparables = NULL) {
            index <- table %in% x
            if (any(index)) {
              which.max(index)
            } else if (missing("nomatch")) {
              NA_integer_
            } else {
              nomatch <- as.integer(nomatch)
              if (length(nomatch) != 1L) {
                NA_integer_
              } else {
                nomatch
              }
            }
          })

#' @rdname match.multiPhylo
#' @aliases match,multiPhylo,multiPhylo-method
#' @export
setMethod("match",
          signature(x = "multiPhylo", table = "multiPhylo"),
          function(x, table, nomatch = NA_integer_, incomparables = NULL) {
            if (missing("nomatch")) {
              nomatch <- NA_integer_
            }
            nomatch <- as.integer(nomatch)
            if (length(nomatch) != 1L) {
              nomatch <- NA_integer_
            }
            
            vapply(x, function(i) {
              match(i, table, nomatch = nomatch)
            }, integer(1))
          })



#' @rdname match.multiPhylo
#' @export
#' @keywords methods
setMethod("%in%",
          signature(x = "multiPhylo", table = "multiPhylo"),
          function(x, table) {
            vapply(x, function (i) any(i %in% table), logical(1))
          })

#' @rdname match.multiPhylo
#' @export
#' @keywords methods
setMethod("%in%",
          signature(x = "multiPhylo", table = "phylo"),
          function(x, table) {
            vapply(x, function (i) all.equal(table, i), logical(1))
          }
          )

#' @rdname match.multiPhylo
#' @export
#' @keywords methods
setMethod("%in%",
          signature(x = "phylo", table = "multiPhylo"),
          function(x, table) {
            any(table %in% x)
          })

#' @rdname match.multiPhylo
#' @export
#' @importFrom ape all.equal.phylo
#' @keywords methods
setMethod("%in%",
          signature(x = "phylo", table = "phylo"),
          function(x, table) all.equal(x, table))

#' @rdname match.Splits
#' @param x,table Splits objects
#' @param return Which index to return: in `x`, in `table`, or both
#' @return `FirstMatchingSplit()` returns an integer
#'  (or length-2 integer if return = "both") specifying the first split in `x`
#'  to have a match in `table`, and the index of that match.
#'  No match is denoted `0` by default.
#' @export
FirstMatchingSplit <- function(x, table, nomatch,
                               return = c("x", "table", "both")) {
  if (!inherits(x, "Splits")) {
    stop("`x` must be a Splits object; try as.Splits(x)")
  }
  if (!inherits(table, "Splits")) {
    stop("`table` must be a Splits object; try as.Splits(table)")
  }
  ij <- first_matching_split_pair(x, table)
  
  if (!missing(nomatch)) {
    ij[ij == 0] <- nomatch
  }
  
  # Return:
  return <- match.arg(return)
  switch(return,
         x     = ij[[1L]],
         table = ij[[2L]],
         both  = ij)
}
