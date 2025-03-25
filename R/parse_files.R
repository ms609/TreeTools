#' Read modification time from "ape" Nexus file
#'
#' `ApeTime()` reads the time that a tree written with "ape" was modified,
#' based on the comment in the Nexus file.
#'
#' @param filepath Character string specifying path to the file.
#' @param format Format in which to return the time: "double" as a sortable numeric;
#'               any other value to return a string in the format
#'               `YYYY-MM-DD hh:mm:ss`.
#'
#' @return `ApeTime()` returns the time that the specified file was created by
#' ape, in the format specified by `format`.
#' @template MRS
#' @export
ApeTime <- function(filepath, format = "double") {
  if (length(filepath) > 1L) {
    stop("`filepath` must be a character string of length 1")
  }
  comment <- readLines(filepath, n = 2)[2]
  Month <- function(month) {
    months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    whichMonth <- months == month
    if (any(whichMonth)) {
      formatC(which(whichMonth), width = 2, flag = "0")
    } else {
      month                                                                     # nocov
    }
  }
  DATEEXP <- ".*? (\\w+)\\s(\\d+)\\s(\\d+\\:\\d\\d\\:\\d\\d)\\s(\\d\\d\\d\\d).*"
  time <- paste0(gsub(DATEEXP, "\\4-", comment),
                 Month(gsub(DATEEXP, "\\1", comment)),
                 gsub(DATEEXP, "-\\2 \\3", comment))

  # Return:
  ifelse(format == "double", as.numeric(as.POSIXct(time, tz = "GMT")), time)
}

#' Extract taxa from a matrix block
#'
#' Extract leaf labels and character states from a Nexus-formatted matrix.
#'
#' @param matrixLines Character vector containing lines of a file that include
#' a phylogenetic matrix. See [`ReadCharacters()`] for expected format.
#' @template characterNumParam
#' @param continuous Logical specifying whether characters are continuous.
#' Treated as discrete if `FALSE`.
#'
#' @return `ExtractTaxa()` returns a matrix with _n_ rows, each named for the
#' relevant taxon, and _c_ columns,
#' each corresponding to the respective character specified in `character_num`.
#'
#' @examples
#' fileName <- paste0(system.file(package = "TreeTools"),
#'                    "/extdata/input/dataset.nex")
#' matrixLines <- readLines(fileName)[6:11]
#' ExtractTaxa(matrixLines)
#'
#' @keywords internal
#' @export
ExtractTaxa <- function(matrixLines, character_num = NULL,
                         continuous = FALSE) {
  taxonLine.pattern <- "('([^']+)'|\"([^\"+])\"|(\\S+))\\s+(.+)$"

  taxonLines <- regexpr(taxonLine.pattern, matrixLines, perl = TRUE) > -1
  # If a line does not start with a taxon name, join it to the preceding line
  taxonLineNumber <- which(taxonLines)
  previousTaxon <- vapply(which(!taxonLines), function(x) {
    max(taxonLineNumber[taxonLineNumber < x])
  }, integer(1))


  taxa <- sub(taxonLine.pattern, "\\2\\3\\4", matrixLines, perl = TRUE)
  taxa <- gsub(" ", "_", taxa, fixed=TRUE)
  taxa[!taxonLines] <- taxa[previousTaxon]
  uniqueTaxa <- unique(taxa)

  tokens <- sub(taxonLine.pattern, "\\5", matrixLines, perl = TRUE)
  if (continuous) {
    tokens <- strsplit(tokens, "\\s+")
    lengths <- lengths(tokens)
    if (length(unique(lengths)) != 1) {
      stop("Different numbers of tokens in different taxa: ",                   # nocov
           paste(lengths, collapse = ", "))                                     # nocov
    }
    tokens <- t(vapply(tokens, I, tokens[[1]]))
  } else {
    tokens <- gsub("\t", "", gsub(" ", "", tokens, fixed = TRUE), fixed = TRUE)
    tokens <- vapply(uniqueTaxa,
                     function(taxon) paste0(tokens[taxa == taxon],
                                             collapse = ""),
                     character(1))
    tokens <- NexusTokens(tokens, character_num = character_num)
  }

  rownames(tokens) <- uniqueTaxa

  # Return:
  tokens
}

#' @rdname ExtractTaxa
#' @param tokens Vector of character strings corresponding to phylogenetic
#'  tokens.
#' @return `NexusTokens()` returns a character vector in which each entry
#' corresponds to the states of a phylogenetic character, or a list containing
#' an error message if input is invalid.
#' @examples
#' NexusTokens("01[01]-?")
#' @export
NexusTokens <- function(tokens, character_num = NULL) {
  tokens.pattern <- "\\([^\\)]+\\)|\\[[^\\]]+\\]|\\{[^\\}]+\\}|\\S"
  matches <- gregexpr(tokens.pattern, tokens, perl = TRUE)

  nChar <- length(matches[[1]])

  if (!exists("character_num") || is.null(character_num)) {
    character_num <- seq_len(nChar)
  } else if (any(character_num > nChar) || any(character_num < 1)) {
    character_num[character_num < 1] <- 1L
    character_num[character_num > nChar] <- nChar
    return(list(paste0("Character number must be between 1 and ", nChar, ".")))
  }

  tokens <- t(vapply(regmatches(tokens, matches),
                     function(x) x[character_num, drop=FALSE],
                     character(length(character_num))))
  if (length(character_num) == 1) {
    tokens <- t(tokens)
  } else if (length(character_num) == 0) {
    stop("No characters selected")
  }

  # Return:
  tokens
}

#' Read phylogenetic characters from file
#'
#' Parse a Nexus \insertCite{Maddison1997}{TreeTools} or 
#' TNT \insertCite{Goloboff2008}{TreeTools} file, reading character states and
#' names.
#'
#' Tested with matrices downloaded from [MorphoBank](https://morphobank.org/)
#' \insertCite{OLeary2011}{TreeTools}, but should also work more widely; please
#' [report](https://github.com/ms609/TreeTools/issues/new?title=Error+parsing+Nexus+file&body=<!--Tell+me+more+and+attach+your+file...-->)
#' incompletely or incorrectly parsed files.
#'
#' Matrices must contain only continuous or only discrete characters;
#' maximum one matrix per file.  Continuous characters will be read as strings
#' (i.e. base type "character").
#'
#' The encoding of an input file will be automatically determined by R.
#' Errors pertaining to an `invalid multibyte string` or
#' `string invalid at that locale` indicate that R has failed to detect
#' the appropriate encoding.  Either
#' [re-save the file](
#' https://support.posit.co/hc/en-us/articles/200532197-Character-Encoding-in-the-RStudio-IDE)
#' in a supported encoding (`UTF-8` is a good choice) or
#' specify the file encoding (which you can find by, for example, opening in
#' [Notepad++](https://notepad-plus-plus.org/downloads/) and identifying
#' the highlighted option in the "Encoding" menu) following the example below.
#'
#' @param filepath character string specifying location of file, or a
#' [connection][base::connections] to the file.
#' @param type Character vector specifying categories of data to extract from
#' file. Setting `type = c("num", "dna")` will return only characters
#' following a `&[num]` or `&[dna]` tag in a TNT input file, listing `num`
#' character blocks before `dna` characters.
#' Leave as `NULL` (the default) to return all characters in their original
#' sequence.
#' @template characterNumParam
#' @param encoding Character encoding of input file.
#'
#' @return `ReadCharacters()` and `ReadTNTCharacters()` return a matrix whose
#' row names correspond to tip labels, and
#' column names correspond to character labels, with the
#' attribute `state.labels` listing the state labels for each character; or
#' a list of length one containing a character string explaining why the
#' function call was unsuccessful.
#'
#' `ReadAsPhyDat()` and `ReadTntAsPhyDat()` return a `phyDat` object.
#'
#' @references \insertAllCited{}
#'
#' @examples
#' fileName <- paste0(system.file(package = "TreeTools"),
#'                    "/extdata/input/dataset.nex")
#' ReadCharacters(fileName)
#'
#' fileName <- paste0(system.file(package = "TreeTools"),
#'                    "/extdata/tests/continuous.nex")
#'
#' continuous <- ReadCharacters(fileName, encoding = "UTF8")
#'
#' # To convert from strings to numbers:
#' at <- attributes(continuous)
#' continuous <- suppressWarnings(as.numeric(continuous))
#' attributes(continuous) <- at
#' continuous
#' @template MRS
#'
#' @seealso
#' - Convert between matrices and `phyDat` objects: [`MatrixToPhyDat()`]
#'
#' - Write characters to TNT-format file: [`WriteTntCharacters()`]
#' @export
ReadCharacters <- function(filepath, character_num = NULL, encoding = "UTF8") {

  lines <- .UTFLines(filepath, encoding)
  nexusComment.pattern <- "\\[[^\\]*?\\]"
  lines <- gsub(nexusComment.pattern, "", lines)
  lines <- trimws(lines)
  lines <- lines[lines != ""]

  semicolons <- which(RightmostCharacter(lines) == ";")
  upperLines <- toupper(lines)

  continuous <- length(grep("DATATYPE[\\S\\=]+CONTINUOUS", upperLines)) > 0
  matrixStart <- which(upperLines == "MATRIX")
  if (length(matrixStart) == 0) {
    return(list("MATRIX block not found in Nexus file."))
  } else if (length (matrixStart) > 1) {
    return(list("Multiple MATRIX blocks found in Nexus file."))                 # nocov
  } else {

    matrixEnd <- semicolons[semicolons > matrixStart][1]
    if (lines[matrixEnd] == ";") matrixEnd <- matrixEnd - 1

    matrixLines <- lines[(matrixStart + 1):matrixEnd]
    tokens <- ExtractTaxa(matrixLines, character_num = character_num,
                          continuous = continuous)

    nChar.pattern <- "DIMENSIONS\\s+.*NCHAR\\s*=\\s*(\\d+)"
    nCharLines <- grepl(nChar.pattern, upperLines, perl = TRUE)
    nChar <- .RegExpMatches("DIMENSIONS\\s+.*NCHAR\\s*=\\s*(\\d+)", upperLines[nCharLines])
    if (length(unique(nChar)) > 1) {
      return(list("Inconsistent DIMENSIONS NCHAR= counts in Nexus file."))
    }
    nChar <- as.integer(nChar[1])
    if (is.null(character_num)) {
      character_num <- seq_len(ncol(tokens))
    }

    ## Written with MorphoBank format in mind: each label on separate line,
    ## each character introduced by integer and terminated with comma.
    labelStart <- grep("^\\s*CHARLABELS\\s*$", upperLines,
                       ignore.case = TRUE, perl = TRUE)
    if (length(labelStart) == 1) {
      labelEnd <- semicolons[semicolons > labelStart][1]
      if (lines[labelEnd] == ";") labelEnd <- labelEnd - 1
      #attr(dat, "char.labels")
      colnames(tokens) <- Unquote(
        .UnescapeQuotes(lines[labelStart + character_num])
      )
    } else {
      if (length(labelStart) > 1)
        return(list("Multiple CharLabels blocks in Nexus file."))
    }

    stateStart <- grep("^\\s*STATELABELS\\s*$", upperLines,
                                 ignore.case = TRUE, perl = TRUE)
    if (length(stateStart) == 1) {
      stateEnd <- semicolons[semicolons > stateStart][1]
      if (is.na(stateEnd)) {
        stop("STATELABELS block missing closing semicolon;")
      }
      stateLines <- lines[stateStart:stateEnd]
      stateStarts <- grep("^\\d+", stateLines)
      stateEnds <- grep("[,;]$", stateLines)
      if (length(stateStarts) != length(stateEnds)) {
        warning("Could not parse character states; does each end with a ' or ;?.")
      } else {
        if (length(character_num) != length(stateStarts)) {
          stateNos <- as.integer(stateLines[stateStarts])
          if (all(!is.na(stateNos))) {
            warning("Missing character state definition for: ",
                    paste0(setdiff(character_num, stateNos), collapse = ", "))
          } else {
            warning("More characters than character state definitions.")
          }
        }

        attr(tokens, "state.labels") <-
          lapply(character_num[seq_along(stateStarts)], function(i)
            Unquote(
              .UnescapeQuotes(
                stateLines[(stateStarts[i] + 1):(stateEnds[i] - 1)]
              )
            )
          )
      }
    } else {
      if (length(labelStart) > 1) {
        return(list("Multiple StateLabels blocks in Nexus file."))
      }
    }

    charStateLabelsStart <- grep("^\\s*CHARSTATELABELS\\s*$", upperLines,
                                 ignore.case = TRUE, perl = TRUE)
    endCommand <- grep(";$", upperLines, perl = TRUE)

    if (length(charStateLabelsStart)) {
      if (length(charStateLabelsStart) > 1) {
        return(list("Multiple CHARSTATELABELS blocks found in Nexus file."))
      }
      charStateLabelsEnd <- endCommand[endCommand > charStateLabelsStart][1]
      if (length(charStateLabelsEnd) == 0) {
        return(list("Unterminated CHARSTATELABELS block in Nexus file."))
      }

      # csl: charStateLabels
      csl <- lines[(charStateLabelsStart + 1):charStateLabelsEnd]
      quote.pattern <- "'[^']*'"
      cslEscapes <- unlist(regmatches(csl, gregexpr(quote.pattern, csl, perl = TRUE)))
      cslEscapes <- substr(cslEscapes, 2L, vapply(cslEscapes, nchar, 1) - 1L)

      cslEscaped <- strsplit(paste(gsub(quote.pattern,
                                        "__TREETOOLS_ESCAPE_SEQUENCE__", csl,
                                        perl = TRUE), collapse = "  "),
                             ",\\s+\\d+", perl = TRUE)[[1]]
      nCsl <- length(cslEscaped)
      if (nCsl != nChar) {
        return(list(paste0("CHARSTATELABELS length (", nCsl,
                           ") does not match NCHAR (", nChar, ")")))
      }
      cslEscaped[1] <- sub("^\\s*1\\b", "", cslEscaped[1], perl = TRUE)
      cslEscaped[nCsl] <- sub("\\s*;\\s*$", "", cslEscaped[nCsl], perl = TRUE)
      cslEscaped <- trimws(sub("\\s*/\\s*", "__TREETOOLS_ESCAPE_SPLITTER__", cslEscaped, perl = TRUE))
      cslEscaped <- gsub("\\s+", "__TREETOOLS_ESCAPE_LABEL_SPLITTER__", cslEscaped, perl = TRUE)
      cslEscaped <- paste(cslEscaped, collapse = "__TREETOOLS_ESCAPE_CHAR_SPLITTER__")
      for (escape in cslEscapes) {
        cslEscaped <- sub("__TREETOOLS_ESCAPE_SEQUENCE__", escape, cslEscaped, fixed = TRUE)
      }
      cslEscaped <- strsplit(cslEscaped, "__TREETOOLS_ESCAPE_CHAR_SPLITTER__", fixed = TRUE)[[1]]
      cslEscaped <- do.call(rbind, strsplit(cslEscaped, "__TREETOOLS_ESCAPE_SPLITTER__", fixed = TRUE))
      colnames(tokens) <- gsub("_", " ", cslEscaped[, 1])[character_num]
      attr(tokens, "state.labels") <- lapply(
        strsplit(cslEscaped[character_num, 2], "__TREETOOLS_ESCAPE_LABEL_SPLITTER__", fixed = TRUE),
        gsub, pattern = "_", replacement = " ", fixed = TRUE)
    }
  }

  # Return:
  tokens
}


#' @rdname ReadCharacters
#' @export
ReadTntCharacters <- function(filepath, character_num = NULL,
                               type = NULL, encoding = "UTF8") {

  lines <- .UTFLines(filepath, encoding)
  tntComment.pattern <- "'[^']*'"
  lines <- gsub(tntComment.pattern, "", lines, perl = TRUE)
  multilineComments <- grep("'", lines, fixed = TRUE)
  nmlc <- length(multilineComments) / 2
  openComment <- multilineComments[seq_len(nmlc) * 2L - 1L]
  closeComment <- multilineComments[seq_len(nmlc) * 2L]
  lines[openComment] <- gsub("'.*", "", lines[openComment])
  lines[closeComment] <- gsub(".*'", "", lines[closeComment])

  lines <- trimws(lines)
  lines <- lines[lines != ""]

  semicolons <- grep(";", lines, fixed = TRUE)
  upperLines <- toupper(lines)

  xread <- grep("^XREAD\\b", lines, ignore.case = TRUE, perl = TRUE)
  if (length(xread) < 1) return(NULL)
  if (length(xread) > 1) {
    message("Multiple character blocks not yet supported;",
            "contact 'TreeTools' maintainer to request.",
            "Returning first block only.")
    xread <- xread[1]
  }

  xreadEnd <- semicolons[semicolons > xread][1]
  if (lines[xreadEnd] == ";") {
    xreadEnd <- xreadEnd - 1L
  }
  xreadLines <- lines[xread:xreadEnd]
  xDimPos <- regexec("'?\\s*(\\d+)\\s*(\\d+)\\s*$", xreadLines)
  xDimLine <- which(vapply(xDimPos, `[`, 1, 1) > -1)[1]
  dimText <- xreadLines[xDimLine]
  dimHit <- xDimPos[[xDimLine]]
  nChar <- as.integer(substr(dimText, dimHit[2], dimHit[2] +
                               attr(dimHit, "match.length")[2] - 1L))
  nTip <- as.integer(substr(dimText, dimHit[3], dimHit[3] +
                              attr(dimHit, "match.length")[3] - 1L))
  matrixLines <- xreadLines[-seq_len(xDimLine)]

  ctypeLines <- grep("^&\\[[\\w\\s]+\\]$", matrixLines, perl = TRUE)
  if (is.null(type)) {
    if (length(ctypeLines)) matrixLines <- matrixLines[-ctypeLines]
  } else {
    types <- matrixLines[ctypeLines]
    blocks <- lapply(paste0("\\b", type, "\\b"), grep, types, ignore.case = TRUE)
    nBlocks <- lengths(blocks)
    if (any(nBlocks == 0)) {
      message("Tags ", paste0(type[nBlocks == 0], collapse = ", "),
      " not found. Ignored: ", types[!seq_along(types) %fin% unlist(blocks)])
    }
    if (all(nBlocks == 0L)) return(NULL)
    blockSpan <- cbind(ctypeLines + 1L, c(ctypeLines[-1] - 1, length(matrixLines)))
    matrixLines <- matrixLines[unlist(
      apply(blockSpan[unlist(blocks), , drop = FALSE],  1,
            function(x) seq.int(from = x[1], to = x[2])))]
  }

  tokens <- ExtractTaxa(matrixLines, character_num)
  if (nrow(tokens) != nTip) {
    warning("Extracted ", nrow(tokens), " taxa, but TNT file specifies ", nTip, # nocov
            ": please check output and report bugs.")                           # nocov
  }
  labelStart <- which(upperLines == "CHARLABELS")
  if (length(labelStart) == 1) {
    labelEnd <- semicolons[semicolons > labelStart][1]
    if (lines[labelEnd] == ";") labelEnd <- labelEnd - 1
    #attr(dat, "char.labels")
    colnames(tokens) <- lines[labelStart + character_num]
  } else {
    if (length(labelStart) > 1)
      return(list("Multiple CharLabels blocks in Nexus file."))
  }


  labelStart <- which(upperLines == "CNAMES")
  if (length(labelStart) == 1) {
    labelEnd <- semicolons[semicolons > labelStart][1]
    if (lines[labelEnd] == ";") labelEnd <- labelEnd - 1
    charLines <- lines[labelStart + character_num]
    charLine.pattern <- "^\\S+\\s\\d+\\s(\\w+)(.*)\\s*;\\s*$"

    # Character labels
    charNames <- gsub(charLine.pattern, "\\1", charLines, perl = TRUE)
    colnames(tokens) <- gsub("_", " ", charNames, fixed = TRUE)

    # State labels
    stateNames <- gsub(charLine.pattern, "\\2", charLines, perl = TRUE)
    attr(tokens, "state.labels") <- lapply(stateNames, function(line) {
      states <- strsplit(trimws(line), "\\s+", perl = TRUE)[[1]]
      trimws(gsub("_", " ", states, fixed = TRUE))
    })
  } else {
    if (length(labelStart) > 1)
      return(list("Multiple cnames entries in TNT file."))
  }

  # Return:
  tokens
}

#' @rdname ReadCharacters
#' @export
ReadTNTCharacters <- ReadTntCharacters

.UTFLines <- function(filepath, encoding) {
  con <- file(filepath, encoding = encoding)
  on.exit(close(con))
  
  tryCatch(
    # Missing EOL might occur in user-generated file, so warning not helpful
    try1 <- enc2utf8(readLines(con, warn = FALSE)),
    warning = function(e) {
      if (substr(e[["message"]], 0, 39) == 
          "invalid input found on input connection") {
        newEnc <- if (toupper(encoding) %in% c("UTF-8", "UTF8")) {
          "latin1"
        } else {
          "UTF8"
        }
        message("Problem reading characters; trying ", newEnc, " file encoding")
        close(con)
        con <- file(filepath, encoding = newEnc)
        enc2utf8(readLines(con, warn = FALSE))
      } else {
        try1
      }
    }
  )
}

.RegExpMatches <- function(pattern, string, i = TRUE, nMatch = 1) {
  res <- regexpr(pattern, string, perl = TRUE, ignore.case = i)
  matches <- regmatches(string, res)
  res[res < 0] <- NA
  res[!is.na(res)] <-
    gsub(pattern, paste0(paste0("\\", seq_len(nMatch)), collapse = ""),
       matches, ignore.case = i, perl = TRUE)
  res
}

#' @rdname ReadCharacters
#' @return `ReadNotes()` returns a list in which each entry corresponds to a
#' single character, and itself contains a list of with two elements:
#'
#' 1. A single character object listing any notes associated with the character
#' 2. A named character vector listing the notes associated with each taxon
#' for that character, named with the names of each note-bearing taxon.
#'
#' @importFrom stats setNames
#' @export
ReadNotes <- function(filepath, encoding = "UTF8") {
  taxon.pattern <- "^\\s+[\"']?([^;]*?)[\"']?\\s*$"
  nTax.pattern <- "DIMENSIONS\\s+.*NTAX\\s*=\\s*(\\d+)"

  lines <- .UTFLines(filepath, encoding)
  upperLines <- toupper(lines)
  trimUpperLines <- trimws(upperLines)

  notesStart <- which(trimUpperLines == "BEGIN NOTES;")
  endBlocks <- which(trimUpperLines %fin% c("END;", "ENDBLOCK;"))
  taxLabels <- which(trimUpperLines == "TAXLABELS")
  semicolons <- which(trimUpperLines == ";")
  nTaxLines <- grepl(nTax.pattern, trimUpperLines, perl = TRUE)


  if (length(notesStart) == 0) {
    return(list("NOTES block not found in Nexus file."))
  } else if (length(notesStart) > 1) {
    return(list("Multiple NOTES blocks found in Nexus file."))
  } else if (length(taxLabels) > 1) {
    return(list("Multiple TAXLABELS found in Nexus file."))
  } else if (!any(nTaxLines)) {
    return(list("No DIMENSIONS NTAX= statment found in Nexus file."))
  } else {
    if (length(taxLabels) == 0) {
      taxa <- names(ReadAsPhyDat(filepath))
    } else {

      nTax <- .RegExpMatches(nTax.pattern, trimUpperLines[nTaxLines])
      if (length(unique(nTax)) > 1) {
        return(list("Inconsistent DIMENSIONS NTAX= counts in Nexus file."))
      }
      nTax <- as.integer(nTax[1])

      taxaEnd <- semicolons[semicolons > taxLabels][1] - 1L
      taxaLines <- lines[(taxLabels + 1):taxaEnd]
      taxon.matches <- grepl(taxon.pattern, taxaLines, perl = TRUE)
      if (sum(taxon.matches) == nTax) {
        taxa <- gsub(taxon.pattern, "\\1", taxaLines[taxon.matches], perl = TRUE)
        taxa <- gsub(" ", "_", taxa, fixed = TRUE)
      } else {
        taxa <- gsub(taxon.pattern, "\\1", taxaLines[taxon.matches], perl = TRUE)
        taxa <- unlist(strsplit(taxa, " "))
        if (length(taxa) != nTax) {
          return(list(paste0("Mismatch: NTAX=", nTax, ", but ", length(taxa),
                             " TAXLABELS found in Nexus file.")))
        }
      }
    }

    notesEnd <- endBlocks[endBlocks > notesStart][1] - 1L
    notesLines <- lines[(notesStart + 1):notesEnd]
    collapsedLines <- paste0(notesLines, collapse = "\r\n")
    # Remove [comments]
    collapsedLines <- gsub("(;\\s*)\\[[^\\]]*\\]", "\\1", collapsedLines,
                           perl = TRUE)
    notes <- strsplit(collapsedLines,
                      # (?i) makes perl regexp case insensitive
                      "(?i)\\r\\n\\s*TEXT\\s+", perl = TRUE)[[1]]

    noteTaxon <- as.integer(.RegExpMatches("\\bTAXON\\s*=\\s*(\\d+)", notes))
    noteChar <- as.integer(.RegExpMatches("\\bCHARACTER\\s*=\\s*(\\d+)", notes))
    noteText <- EndSentence(MorphoBankDecode(
      .RegExpMatches("\\bTEXT\\s*=\\s*['\"]([\\s\\S]+)['\"];\\s*$", notes)))

    seqAlongNotes <- if (any(!is.na(noteChar))) {
      seq_len(max(noteChar, na.rm = TRUE))
    } else {
      numeric(0)
    }

    # Return:
    setNames(lapply(seqAlongNotes, function(i) {
      byTaxon <- !is.na(noteChar) & noteChar == i & !is.na(noteTaxon)
      ret <- list(
        noteText[!is.na(noteChar) & noteChar == i & is.na(noteTaxon)],
        setNames(noteText[byTaxon], taxa[noteTaxon[byTaxon]]))

      # Return:
      ret
    }), seqAlongNotes)
  }
}

#' Add full stop to end of a sentence
#'
#' @param string Input string
#'
#' @return `EndSentence()` returns `string`, punctuated with a final full stop
#' (period).`
#'
#' @examples
#' EndSentence("Hello World") # "Hello World."
#' @author Martin R. Smith
#' @family string parsing functions
#' @export
EndSentence <- function(string) {
  if (length(string)) {
    ret <- gsub("\\s*\\.?\\s*\\.$", ".", paste0(string, "."), perl = TRUE)
    ret <- gsub("(\\.[\"'])\\.$", "\\1", ret, perl = TRUE)
    ret <- gsub("([!\\?])\\.$", "\\1", ret, perl = TRUE)
    ret[ret == "."] <- ""
    ret
  } else {
    string
  }
}

#' Remove quotation marks from a string
#'
#' @param string Input string
#'
#' @return `Unquote()` returns `string`, with any matched punctuation marks
#' and trailing whitespace removed.
#'
#' @examples
#' Unquote("'Hello World'")
#' @author Martin R. Smith
#' @family string parsing functions
#' @export
Unquote <- function(string) {
  noSingle <- vapply(string, gsub, character(1),
                     pattern = "^\\s*'\\s*(.*?)\\s*'\\s*$",
                     replacement = "\\1", USE.NAMES = FALSE)
  vapply(noSingle, gsub, character(1),
         pattern = "^\\s*\"\\s*(.*?)\\s*\"\\s*$", replacement = "\\1",
         USE.NAMES = FALSE)
}

# Unescape quotes in Nexus format
.UnescapeQuotes <- function(string) {
  gsub("(?=.)''(?=.)", "'", string, perl = TRUE)
}

#' Decode MorphoBank text
#'
#' Converts strings from MorphoBank notes into a Latex-compatible format.
#'
#' @param string String to process
#'
#' @return `MorphoBankDecode()` returns a string with new lines and punctuation
#' reformatted.
#' @family string parsing functions
#' @author Martin R. Smith
#' @export
MorphoBankDecode <- function(string) {
  string <- gsub("^n", "  \n", string, fixed = TRUE)
  string <- gsub("''", "'", string, fixed = TRUE)
  string <- gsub(" - ", " -- ", string, fixed = TRUE)
  string <- gsub("(\\d)\\-(\\d)", "\\1--\\2", string, perl = TRUE)
  string <- gsub("(\\d) ?um\\b", "\\1 \u{03BC}m", string, perl = TRUE)
  string <- gsub(" [recoded as neomorphic]",
                 " Inapplicable tokens in this neomorphic character have been replaced with the absent token, following @Brazeau2018", string, fixed = TRUE)

  # Return:
  string
}

#' Convert between matrices and `phyDat` objects
#'
#' `MatrixToPhyDat()` converts a matrix of tokens to a `phyDat` object;
#' `PhyDatToMatrix()` converts a `phyDat` object to a matrix of tokens.
#'
#' @param tokens Matrix of tokens, possibly created with [`ReadCharacters()`]
#' or [`ReadTntCharacters()`].
#' Row names should correspond to leaf labels; column names may optionally
#' correspond to character labels.
#'
#' @return `MatrixToPhyDat()` returns an object of class `phyDat`.
#'
#' @family phylogenetic matrix conversion functions
#' @examples 
#' tokens <- matrix(c(0, 0, "0", 0, 0,
#'                    0, 0, "1", 0, 1,
#'                    0, 0, "1", 0, 1,
#'                    0, 0, "2", 0, 1,
#'                    1, 1, "-", 1, 0,
#'                    1, 1, "2", 1, "{01}"),
#'                    nrow = 6, ncol = 5, byrow = TRUE,
#'                    dimnames = list(
#'                      paste0("Taxon_", LETTERS[1:6]),
#'                      paste0("Char_", 1:5)))
#'                    
#' MatrixToPhyDat(tokens)
#' @template MRS
#' @export
MatrixToPhyDat <- function(tokens) {
  if (inherits(tokens, "phyDat")) {
    warning("MatrixToPhyDat() expects a matrix, not a phyDat object")
    return(tokens)
  }
  allTokens <- unique(as.character(tokens))
  if (any(nchar(allTokens) == 0)) {
    problems <- apply(tokens, 1, function(x) which(nchar(x) == 0))
    problemTaxa <- lengths(problems) > 0
    problemTaxa <- names(problemTaxa[problemTaxa])
    warning("Blank tokens ('') found in taxa: ",
            paste0(problemTaxa, collapse = ", "))
  }
  tokenNumbers <- seq_along(allTokens)
  names(tokenNumbers) <- allTokens
  matches <- gregexpr("[\\d\\-\\w]", allTokens, perl = TRUE)
  whichTokens <- regmatches(allTokens, matches)
  levels <- sort(unique(unlist(whichTokens)))
  whichTokens[allTokens == "?"] <- list(levels)
  contrast <- vapply(whichTokens, function(x) levels %fin% x,
                     logical(length(levels)))
  contrast <- 1 * if (is.null(dim(contrast))) {
    as.matrix(contrast)
  } else {
    t(contrast)
  }
  dimnames(contrast) <- list(allTokens, levels)
  dat <- .PhyDatWithContrast(tokens, contrast = contrast)

  # Return:
  dat
}


#' @export
`[.phyDat` <- .SubsetPhyDat <- function(x, i, j, ..., drop = FALSE) {
  mat <- PhyDatToMatrix(x)
  MatrixToPhyDat(mat[i, j, ..., drop = FALSE])
}


.PhyDatWithContrast <- function(dat, contrast) {
  if (is.null(dim(dat))) {
    dat <- t(t(dat))
  }
  tipLabels <- rownames(dat)
  if (is.null(tipLabels)) {
    stop("Data rows must be named with tip labels")
  }
  
  levelSet <- dimnames(contrast)
  allLevels <- levelSet[[1]]
  levels <- levelSet[[2]]
  rownames(contrast) <- NULL
  
  # See https://stackoverflow.com/questions/70557817
  groups <- do.call(grouping, as.data.frame(t(dat)))
  ends <- attr(groups, "ends")
  i <- rep(seq_along(ends), c(ends[1], diff(ends)))[order(groups)]
  firstOccurrence <- match(i, i)
  tab <- table(firstOccurrence)
  weight <- as.integer(tab)
  tab[] <- seq_along(tab)
  index <- as.integer(tab[as.character(firstOccurrence)])
  
  duplicate <- duplicated(firstOccurrence)
  phyMat <- matrix(match(dat[, !duplicate], allLevels),
                   dim(dat)[1], sum(!duplicate))
  
  # Return:
  structure(
    lapply(seq_along(tipLabels), function(i) phyMat[i, ]),
    #as.list(asplit(phyMat, 1)),
    names  = tipLabels,
    weight = as.integer(weight),
    nr = length(weight),
    nc = length(levels),
    index = unname(index),
    levels = levels,
    allLevels = allLevels,
    type = "USER",
    contrast = contrast,
    class = "phyDat")
}


#' @rdname MatrixToPhyDat
#' @param dataset A dataset of class `phyDat`.
#' @param ambigNA,inappNA Logical specifying whether to denote ambiguous /
#' inapplicable characters as `NA` values.
#' @param parentheses Character vector specifying style of parentheses
#' with which to enclose ambiguous characters. `c("[", "]")` or `"[]"` will
#' render `[01]`.
#' `NULL` will use the token specified in the `phyDat` object; but beware that
#' this will be treated as a distinct (non-ambiguous) token if re-encoding with
#' `PhyDatToMatrix()`.
#' @param sep Character with which to separate ambiguous tokens, e.g. `','`
#' will render `[0,1]`.
#' @return `PhyDatToMatrix()` returns a matrix corresponding to the
#' uncompressed character states within a `phyDat` object.
#' @examples 
#' data("Lobo", package = "TreeTools")
#' head(PhyDatToMatrix(Lobo.phy)[, 91:93])
#' @export
PhyDatToMatrix <- function(dataset, ambigNA = FALSE, inappNA = ambigNA,
                           parentheses = c("{", "}"), sep = "") {
  if (!is.null(parentheses) && length(parentheses) == 0) {
    parentheses <- c("", "")
  } else if (length(parentheses) == 1) {
    parentheses <- c(
      substr(parentheses, 1, 1),
      substr(parentheses, nchar(parentheses), nchar(parentheses))
    )
  }
  
  at <- attributes(dataset)
  allLevels <- as.character(at[["allLevels"]])
  if (inappNA) {
    allLevels[allLevels == "-"] <- NA_character_
  }
  if (ambigNA) {
    allLevels[rowSums(at[["contrast"]]) != 1L] <- NA_character_
  } else if (!is.null(parentheses)) {
    cont <- at[["contrast"]]
    nTokens <- rowSums(cont)
    levels <- colnames(cont)
    partAmbig <- nTokens != 1L & nTokens < dim(cont)[2]
    allLevels[partAmbig] <- paste0(
      parentheses[1],
      apply(cont[partAmbig, , drop = FALSE] > 0, 1, function(x) {
        paste0(levels[x], collapse = sep)
      }),
      parentheses[2])
  }
  matrix(allLevels[unlist(dataset, recursive = FALSE, use.names = FALSE)],
         ncol = at[["nr"]], byrow = TRUE, dimnames = list(at[["names"]], NULL)
         )[, at[["index"]], drop = FALSE]
}

#' @rdname ReadCharacters
#' @export
ReadAsPhyDat <- function(...) {
  MatrixToPhyDat(ReadCharacters(...))
}


#' @rdname ReadCharacters
#' @param \dots Parameters to pass to `Read[Tnt]Characters()`.
#' @export
ReadTntAsPhyDat <- function(...) {
  MatrixToPhyDat(ReadTntCharacters(...))
}

#' @rdname ReadCharacters
#' @export
ReadTNTAsPhyDat <- ReadTntAsPhyDat


#' @describeIn ReadCharacters A convenient wrapper for \pkg{phangorn}'s
#' `phyDat()`, which converts a **list** of morphological characters into a
#' `phyDat` object.
#' If your morphological characters are in the form of a **matrix**, perhaps
#' because they have been read using [`read.table()`], try [`MatrixToPhyDat()`]
#' instead.
#'
#' @param dataset list of taxa and characters, in the format produced by 
#' \code{\link[ape]{read.nexus.data}()}:
#'   a list of sequences each made of a single character vector,
#'   and named with the taxon name.
#'
#' @export
PhyDat <- function(dataset) {
  nChar <- length(dataset[[1]])
  if (nChar == 1) {
    mat <- matrix(unlist(dataset), dimnames = list(names(dataset), NULL))
  } else {
    mat <- t(vapply(dataset, I, dataset[[1]]))
  }
  MatrixToPhyDat(mat)
}

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
#' the string.
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
  ret <- vapply(phy,
                function(x) levelTranslation[x[phyIndex]],
                character(length(phyIndex)))
  ret <- if (concatenate || is.null(dim(ret))) { # If only one row, don't need to apply
    if (!byTaxon) ret <- t(ret)
    paste0(c(ret, ps), collapse = "")
  } else {
    if (byTaxon) ret <- t(ret)
    paste0(apply(ret, 1, paste0, collapse = ""), ps)
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


#' Rightmost character of string
#'
#' `RightmostCharacter()` is a convenience function that returns the final
#' character of a string.
#'
#' @param string Character string.
#' @param len (Optional) Integer specifying number of characters in `string`.
#'
#' @return `RightmostCharacter()` returns the rightmost character of a string.
#' @examples
#' RightmostCharacter("Hello, World!")
#'
#' @template MRS
#' @export
#' @family string parsing functions
RightmostCharacter <- function(string, len = nchar(string)) {
  substr(string, len, len)
}

#' Write Newick Tree
#'
#' `NewickTree()` encodes a tree as a Newick-format string.
#' This differs from [`write.tree()`][ape::write.tree] in the encoding of
#' spaces as spaces, rather than underscores.
#'
#' @template treeParam
#'
#' @return `NewickTree()` returns a character string denoting `tree` in Newick
#' format.
#'
#' @examples
#' NewickTree(BalancedTree(LETTERS[4:9]))
#'
#' @seealso Use tip numbers, rather than leaf labels: [`as.Newick`]
#' @importFrom ape write.tree
#' @export
NewickTree <- function(tree) gsub("_", " ", write.tree(tree), fixed = TRUE)
