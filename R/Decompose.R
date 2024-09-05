#' Decompose additive (ordered) phylogenetic characters
#' 
#' `Decompose()` decomposes additive characters into a series of binary
#' characters, which is mathematically equivalent when analysed under
#' equal weights parsimony.  (This equivalence is not exact
#' under implied weights or under probabilistic tree inference methods.)
#' 
#' An ordered (additive) character can be rewritten as a mathematically
#' equivalent hierarchy of binary neomorphic characters
#' \insertCite{Farris1970}{TreeTools}.  
#' Two reasons to prefer the latter approach are:
#' 
#' - It makes explicit the evolutionary assumptions underlying an ordered 
#'   character, whether the underlying ordering is linear, reticulate or
#'   branched \insertCite{Mabee1989}{TreeTools}.
#' - It avoids having to identify characters requiring special treatment to
#'   phylogenetic software, which requires the maintenance of an up-to-date 
#'   log of which characters are treated as additive and which sequence their
#'   states occur in, a step that may be overlooked by re-users of the data.
#'   
#' Careful consideration is warranted when evaluating whether a group of
#' related characteristics ought to be treated as ordered
#' \insertCite{Wilkinson1992}{TreeTools}.
#' On the one hand, the 'principle of indifference' states that we should treat
#' all transformations as equally probable (/ surprising / informative);
#' ordered characters fail this test, as larger changes are treated as less
#' probable than smaller ones.
#' On the other hand, ordered characters allow more opportunities for homology
#' of different character states, and might thus be defended under the auspices
#' of Hennig’s Auxiliary Principle \insertCite{Wilkinson1992}{TreeTools}.
#' 
#' For a case study of how ordering phylogenetic characters can affect 
#' phylogenetic outcomes in practice, see
#' \insertCite{Brady2024;textual}{TreeTools}.
#'
#' @references \insertAllCited{}
#' @template datasetParam
#' @param indices Integer or logical vector specifying indices of characters
#' that should be decomposed
#' 
#' @return `Decompose()` returns a `phyDat` object in which the specified
#' ordered characters have been decomposed into binary characters.  
#' The attribute `originalIndex` lists the index of the character in
#' `dataset` to which each element corresponds.
#' 
#' @examples
#' data("Lobo")
#' 
#' # Identify character 11 as additive
#' # Character 11 will be replaced with two characters
#' # The present codings 0, 1 and 2 will be replaced with 00, 10, and 11.
#' decomposed <- Decompose(Lobo.phy, 11)
#' 
#' NumberOfChars <- function(x) sum(attr(x, "weight"))
#' NumberOfChars(Lobo.phy)   # 115 characters in original
#' NumberOfChars(decomposed) # 116 characters in decomposed
#' @template MRS
#' @family phylogenetic matrix conversion functions
#' @export
Decompose <- function(dataset, indices) {
  if (is.matrix(dataset)) {
    dataset <- MatrixToPhyDat(dataset)
  }
  if (!inherits(dataset, "phyDat")) {
    stop("`dataset` must be an object of class `phyDat`")
  }
  levels <- attr(dataset, "levels")
  inappLevel <- levels == "-"
  appLevels <- levels[!inappLevel]
  mat <- as.matrix(dataset)
  nTaxa <- dim(mat)[[1]]
  nChar <- dim(mat)[[2]]
  if (is.numeric(indices)) {
    if (any(indices > nChar)) {
      warning("Only ", nChar, " characters; indices ",
              paste(indices[indices > nChar], collapse = ", "),
              " not found.")
    }
    indices <- seq_len(nChar) %in% indices
  } else {
    if (length(indices) != nChar) {
      stop("`length(indices)` must correspond to the number of characters in ",
           "`dataset`")
    }
  }
  
  replacements <- if (getRversion() >= "4.1") {
    # apply(simplify = FALSE) becomes available in R4.1.0
    apply(mat[, indices, drop = FALSE], 2, function(char) {
      whichLevels <- which(vapply(appLevels,
                                 function(x) any(grepl(x, char, fixed = TRUE)),
                                 logical(1)))
      if (!length(whichLevels)) {
        return(matrix(character(0), nrow = nTaxa, ncol = 0))
      }
      maxLevel <- max(whichLevels)
      vapply(seq_len(maxLevel)[-1], function(i) {
        zero <- .RegExpEscape(appLevels[1:(i - 1)])
        one <- .RegExpEscape(appLevels[(i):maxLevel])
        gsub("(0)0+|(1)1+", "\\1",
             gsub(paste0(c("[", one, "]"), collapse = ""), "1",
                  gsub(paste0(c("[", zero, "]"), collapse = ""), "0", char)
             )
        )
      }, char)
    }, simplify = FALSE)
  } else {                                                         # nocov start
    lapply(which(indices), function(i) {
      char <- mat[, i]
      whichLevels <- which(vapply(appLevels,
                                 function(x) any(grepl(x, char, fixed = TRUE)),
                                 logical(1)))
      if (!length(whichLevels)) {
        return(matrix(character(0), nrow = nTaxa, ncol = 0))
      }
      maxLevel <- max(whichLevels)
      vapply(seq_len(maxLevel)[-1], function(i) {
        zero <- .RegExpEscape(appLevels[1:(i - 1)])
        one <- .RegExpEscape(appLevels[(i):maxLevel])
        gsub("(0)0+|(1)1+", "\\1",
             gsub(paste0(c("[", one, "]"), collapse = ""), "1",
                  gsub(paste0(c("[", zero, "]"), collapse = ""), "0", char)
             )
        )
      }, char)
    })
  }                                                                # nocov end
  nNew <- `[<-`(rep(1, nChar), indices, vapply(replacements, ncol, 1))
  ret <- matrix(nrow = nTaxa, ncol = sum(nNew),
                dimnames = list(dimnames(mat)[[1]], NULL))
  newInd <- cumsum(nNew)
  ret[, newInd[!indices]] <- mat[, !indices]
  
  offset <- c(0, cumsum(nNew[indices] - 1))
  
  whichInds <- which(indices)
  replInd <- unlist(lapply(seq_along(whichInds), function(i) {
    index <- whichInds[[i]]
    index + offset[[i]] + seq_len(nNew[[index]]) - 1
    }
  ))
  if (!is.null(replInd)) {
    ret[, replInd] <- do.call(cbind, replacements)
  }
  
  # Consistent ambiguity symbols
  ret <- gsub("\\[(\\d)\\]", "\\1", 
              gsub("{", "[", fixed = TRUE,
                   gsub("}", "]", fixed = TRUE, ret)
              )
            )
  
  # Return:
  structure(
    MatrixToPhyDat(ret),
    originalIndex = rep(seq_len(nChar), nNew)
  )
}

.RegExpEscape <- function(x) {
  gsub("([\\[\\]\\{\\}\\?\\.\\+\\*\\|\\^\\$\\-\\(\\)\\\\/])", "\\\\\\1", x,
       perl = TRUE)
}
