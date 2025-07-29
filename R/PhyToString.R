#' @rdname PhyToString
#'
#' @param string String of tokens, optionally containing whitespace, with no
#'   terminating semi-colon.
#' @param tips (Optional) Character vector corresponding to the names (in order)
#' of each taxon in the matrix, or an object such as a tree from which
#' tip labels can be extracted.
#' @param byTaxon Logical; if `TRUE`, string is one **taxon's** coding at a
#' time; if `FALSE`, string is interpreted as one **character's** coding at a
#' time.
#'
#' @return `StringToPhyDat()` returns an object of class `phyDat`.
#'
#' @examples
#' StringToPhyDat("-?01231230?-", c("Lion", "Gazelle"), byTaxon = TRUE)
#' # encodes the following matrix:
#' # Lion     -?0123
#' # Gazelle  1230?-
#'
#' @export
StringToPhyDat <- function(string, tips, byTaxon = TRUE) {
  tokens <- NexusTokens(string)
  if (missing(tips)) {
    tips <- length(tokens)
  }
  tips <- TipLabels(tips)
  tokens <- matrix(tokens, nrow = length(tips), byrow = byTaxon,
                   dimnames = list(tips, NULL))
  
  # Return:
  MatrixToPhyDat(tokens)
}
#' @rdname PhyToString
StringToPhydat <- StringToPhyDat

#' Convert between strings and `phyDat` objects
#'
#' `PhyDatToString()` converts a `phyDat` object as a string;
#' `StringToPhyDat()` converts a string of character data to a `phyDat` object.
#'
#' @param phy An object of class `phyDat`.
#' @param parentheses Character specifying format of parentheses with which to
#' surround ambiguous tokens.  Choose from: \code{\{} (default), `[`, `(`, `<`.
#' @param collapse Character specifying text, perhaps `,`, with which to
#' separate multiple tokens within parentheses.
#' @param ps Character specifying text, perhaps `;`, to append to the end of
#' each string.
#' @param useIndex Logical (default: `TRUE`) specifying whether to print
#' duplicate characters multiple times, as they appeared in the original matrix.
#' @param byTaxon Logical. If `TRUE`, write one taxon followed by the next.
#' If `FALSE`, write one character followed by the next.
#' @param concatenate Logical specifying whether to concatenate all
#' characters/taxa into a single string, or to return a separate string for
#' each entry.
#'
#' @examples
#' fileName <- paste0(system.file(package = "TreeTools"),
#'                    "/extdata/input/dataset.nex")
#' phyDat <- ReadAsPhyDat(fileName)
#' PhyToString(phyDat, concatenate = FALSE)
#'
#' @return `PhyToString()` returns a character vector listing a text
#' representation of the phylogenetic character state for each taxon in turn.
#'
#' @family phylogenetic matrix conversion functions
#' @template MRS
#' @importFrom stringi stri_paste
#' @export
PhyToString <- function(phy, parentheses = "{", collapse = "", ps = "",
                        useIndex = TRUE, byTaxon = TRUE, concatenate = TRUE) {
  at <- attributes(phy)
  phyLevels <- at[["allLevels"]]
  if (sum(phyLevels == "-") > 1) {
    stop("More than one inapplicable level identified.  Is phy$levels malformed?")
  }
  phyChars <- at[["nr"]]
  phyContrast <- at[["contrast"]] == 1
  phyIndex <- if (useIndex) {
    at[["index"]]
  } else {
    seq_len(phyChars)
  }
  outLevels <- at[["levels"]]
  
  levelLengths <- vapply(outLevels, nchar, integer(1))
  longLevels <- levelLengths > 1
  if (any(longLevels)) {
    if ("10" %fin% outLevels && !(0 %fin% outLevels)) {
      outLevels[outLevels == "10"] <- "0"
      longLevels["10"] <- FALSE
    }
    outLevels[longLevels] <- LETTERS[seq_len(sum(longLevels))]
  }
  
  switch(parentheses,
         "(" = {openBracket <- "("; closeBracket = ")"},
         ")" = {openBracket <- "("; closeBracket = ")"},
         "<" = {openBracket <- "<"; closeBracket = ">"},
         ">" = {openBracket <- "<"; closeBracket = ">"},
         "[" = {openBracket <- "["; closeBracket = "]"},
         "]" = {openBracket <- "["; closeBracket = "]"},
         {openBracket <- "{"; closeBracket = "}"})
  
  levelTranslation <- apply(phyContrast, 1, function(x)
    ifelse(sum(x) == 1, as.character(outLevels[x]),
           paste0(c(openBracket, paste0(outLevels[x], collapse = collapse),
                    closeBracket), collapse = ""))
  )
  if (any(ambigToken <- apply(phyContrast, 1, all))) {
    levelTranslation[ambigToken] <- "?"
  }
  ret <- vapply(phy, USE.NAMES = FALSE,
                function(x) levelTranslation[x[phyIndex]],
                character(length(phyIndex)))
  
  ret <- if (isTRUE(concatenate) || is.null(dim(ret))) 
    { # If only one row, don't need to apply
    if (!isTRUE(byTaxon)) {
      ret <- t(ret)
    }
    if (isTRUE(concatenate)) {
      stri_paste(c(ret, ps), collapse = "")
    } else {
      if (isTRUE(byTaxon)) {
        stri_paste(stri_paste(ret, ps))
      } else {
        stri_paste(ret, ps, collapse = "")
      }
    }
  } else {
    if (isTRUE(byTaxon)) ret <- t(ret)
    stri_paste(apply(ret, 1, stri_paste, collapse = ""), ps)
  }
  # Return:
  ret
}

#' @rdname PhyToString
#' @export
PhyDatToString <- PhyToString
#' @rdname PhyToString
#' @export
PhydatToString <- PhyToString
