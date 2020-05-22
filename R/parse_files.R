#' Read modification time from ape Nexus file
#'
#' Reads the time that an ape tree was modified from the comment in the Nexus
#' file.
#'
#' @param filename Character string specifying path to the file
#' @param format Format in which to return the time: 'double' as a sortable numeric;
#'               any other value to return a string in the format YYYY-MM-DD hh:mm:ss
#'
#' @return `ApeTime` returns the time that the specified file was created by
#' ape, in the format specified by `format`.
#' @export
#' @template MRS
#'
ApeTime <- function (filename, format='double') {
  if (length(filename) > 1L) stop("`filename` must be a character string of length 1")
  comment <- readLines(filename, n=2)[2]
  Month <- function (month) {
    months <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
    whichMonth <- months == month
    if (any(whichMonth)) {
      formatC(which (whichMonth), width=2, flag="0")
    } else {
      month
    }
  }
  DATEEXP <- ".*? (\\w+)\\s(\\d+)\\s(\\d+\\:\\d\\d\\:\\d\\d)\\s(\\d\\d\\d\\d).*"
  time <- paste0(gsub(DATEEXP, "\\4-", comment),
                 Month(gsub(DATEEXP, "\\1", comment)),
                 gsub(DATEEXP, "-\\2 \\3", comment))

  # Return:
  ifelse(format=='double', as.numeric(as.POSIXct(time, tz = "GMT")), time)
}

#' Parse TNT Tree
#'
#' Reads a tree from TNT's parenthetical output.
#'
#' [TNT](http://www.lillo.org.ar/phylogeny/tnt/) is software for parsimony
#' analysis.  Whilst its implementation of tree search is extremely rapid,
#' analysis of results in TNT is made difficult by its esoteric and scantly
#' documented scripting language.
#'
#' `ReadTNTTree` aims to aid the user by facilitating the import of trees
#' generated in TNT into R for further analysis.
#'
#' The function depends on tree files being saved by TNT in parenthetical
#' notation, using the TNT command `tsav*`.
#' Trees are easiest to load into R if taxa have been saved using their names
#' (TNT command `taxname=`).  In this case, the TNT `.tre` file
#' contains tip labels and can be parsed directly.  The downside is that the
#' uncompressed `.tre` files will have a larger file size.
#'
#' `ReadTNTTree` can also read `.tre` files in which taxa have been saved using
#' their numbers (`taxname-`).  Such files contain a hard-coded link to the
#' matrix file that was used to generate the trees, in the first line of the
#' `.tre` file.  This poses problems for portability: if the matrix file is
#' moved, or the `.tre` file is accessed on another computer, the taxon names
#' may be lost.  As such, it is important to check that the matrix file
#' exists in the expected location -- if it does not,
#' either use the `relativePath` argument to point to its new location, or
#' specify `tipLabels` to manually specify the tip labels.
#'
#' @param filename character string specifying path to TNT `.tre` file,
#' relative to the R working directory (visible with `getwd()`).
#' @param relativePath (discouraged) character string specifying location of the
#' matrix file used to generate the TNT results, relative to the current working
#' directory.  Taxon names will be read from this file if they are not specified
#'  by `tipLabels`.
#' @param keepEnd (optional, default 1) integer specifying how many elements of
#'  the file path to conserve when creating relative path (see examples).
#' @param tipLabels (optional) character vector specifying the names of the
#' taxa, in the sequence that they appear in the TNT file.  If not specified,
#' taxon names will be loaded from the data file linked in the first line of the
#'  `.tre` file specified in `filename`.
#'
#' @return `ReadTNTTree` returns a tree of class \code{phylo}, corresponding
#' to the tree in `filename`.
#'
#' @examples
#'   # In the examples below, TNT has read a matrix from
#'   # "c:/TreeTools/input/dataset.nex"
#'   # The results of an analysis were written to
#'   # "c:/TreeTools/output/results1.tnt"
#'   #
#'   # results1.tnt will contain a hard-coded reference to
#'   # "c:/TreeTools/input/dataset.nex".
#'
#'   # On the original machine (but not elsewhere), it would be possible to read
#'   # this hard-coded reference from results.tnt:
#'   # ReadTntTree('output/results1.tnt')
#'
#'   # These datasets are provided with the `TreeTools` package, which will
#'   # probably not be located at c:/TreeTools on your machine:
#'
#'   oldWD <- getwd() # Remember the current working directory
#'   setwd(system.file(package = 'TreeTools'))
#'
#'   # If taxon names were saved within the file (using `taxname=` in TNT),
#'   # then our job is easy:
#'   ReadTntTree('extdata/output/named.tre')
#'
#'   # But if taxa were compressed to numbers (using `taxname-`), we need to
#'   # look up the original matrix in order to dereference the tip names.
#'   #
#'   # We need to extract the relevant file path from the end of the
#'   # hard-coded path in the original file.
#'   #
#'   # We are interested in the last two elements of
#'   # c:/TreeTools/input/dataset.nex
#'   #                2      1
#'   #
#'   # '.' means "relative to the current directory"
#'   ReadTntTree('extdata/output/numbered.tre', './extdata', 2)
#'
#'   # If working in a lower subdirectory
#'   setwd('./extdata/otherfolder')
#'
#'   # then it will be necessary to navigate up the directory path with '..':
#'   ReadTntTree('../output/numbered.tre', '..', 2)
#'
#'
#'   setwd(oldWD) # Restore original working directory
#'
#' @template MRS
#' @importFrom ape read.tree
#' @export
ReadTntTree <- function (filename, relativePath = NULL, keepEnd = 1L,
                         tipLabels = NULL) {
  fileText <- readLines(filename)
  trees <- lapply(fileText[2:(length(fileText)-1)], TNTText2Tree)

  if (!any(grepl('[A-z]', trees[[1]]$tip.label))) {
    if (is.null(tipLabels)) {
      taxonFile <- gsub("tread 'tree(s) from TNT, for data in ", '',
                        fileText[1], fixed=TRUE)
      taxonFile <- gsub("'", '', gsub('\\', '/', taxonFile, fixed=TRUE),
                        fixed=TRUE)
      if (!is.null(relativePath)) {
        taxonFileParts <- strsplit(taxonFile, '/')[[1]]
        nParts <- length(taxonFileParts)
        if (nParts < keepEnd) {
          stop("Taxon file path (", taxonFile, ") contains fewer than keepEnd (",
               keepEnd, ") components.")
        }
        taxonFile <- paste0(c(relativePath,
                              taxonFileParts[(nParts + 1L - keepEnd):nParts]),
                            collapse='/')
      }

      if (!file.exists(taxonFile)) {
        warning("Cannot find linked data file:\n  ", taxonFile)
      } else {
        tipLabels <- rownames(ReadTntCharacters(taxonFile, 1))
        if (is.null(tipLabels)) {
          # TNT character read failed.  Perhaps taxonFile is in NEXUS format?
          tipLabels <- rownames(ReadCharacters(taxonFile, 1))
        }
        if (is.null(tipLabels)) {
          warning("Could not read taxon names from linked TNT file:\n  ",
                  taxonFile,
                  "\nIs the file in TNT or Nexus format?",
                  " If failing inexplicably, please report:",
                  "\n  https://github.com/ms609/TreeTools/issues/new")
        }
      }
    }

    trees <- lapply(trees, function (tree) {
      tree$tip.label <- tipLabels[as.integer(tree$tip.label) + 1]
      tree
    })
  }

  # Return:
  if (length(trees) == 1) {
    trees[[1]]
  } else if (length(trees) == 0) {
    NULL
  } else {
    class(trees) <- 'multiPhylo'
    trees
  }

}

#' @describeIn ReadTntTree Converts text representation of a tree in TNT to an
#'  object of class `phylo`.
#' @param treeText Character string describing a tree, in the parenthetical
#'                 format output by TNT.
#' @template MRS
#' @export
TNTText2Tree <- function (treeText) {
  treeText <- gsub("([\\w']+)", "\\1,", treeText, perl = TRUE)
  treeText <- gsub(")(", "),(", treeText, fixed = TRUE)
  treeText <- gsub("*", ";", treeText, fixed = TRUE)
  # Return:
  read.tree(text=gsub(", )", ")", treeText, fixed=TRUE))
}

#' Extract taxa from a matrix block
#'
#' Reads the character information from a Nexus-formatted matrix into R.
#'
#' Extracts the names of tips from the matrix block of a Nexus file.
#'
#' @param matrixLines lines of a file containing a phylogenetic matrix
#'  (see [`ReadCharacters`] for expected format).
#' @template characterNumParam
#' @template sessionParam
#' @param continuous Logical specifying whether characters are continuous.
#' Treated as discrete if `FALSE`.
#'
#' @return `ExtractTaxa` returns a matrix with _n_ rows, each named for the
#' relevant taxon, and _c_ columns,
#' each corresponding to the respective character specified in `character_num`.
#'
#' @examples
#' fileName <- paste0(system.file(package='TreeTools'),
#'                    '/extdata/input/dataset.nex')
#' matrixLines <- readLines(fileName)[6:11]
#' ExtractTaxa(matrixLines)
#'
#' @keywords internal
#' @export
ExtractTaxa <- function (matrixLines, character_num = NULL, session = NULL,
                         continuous = FALSE) {
  taxonLine.pattern <- "('([^']+)'|\"([^\"+])\"|(\\S+))\\s+(.+)$"

  taxonLines <- regexpr(taxonLine.pattern, matrixLines, perl = TRUE) > -1
  # If a line does not start with a taxon name, join it to the preceding line
  taxonLineNumber <- which(taxonLines)
  previousTaxon <- vapply(which(!taxonLines), function (x) {
    max(taxonLineNumber[taxonLineNumber < x])
  }, integer(1))


  taxa <- sub(taxonLine.pattern, "\\2\\3\\4", matrixLines, perl = TRUE)
  taxa <- gsub(" ", "_", taxa, fixed=TRUE)
  taxa[!taxonLines] <- taxa[previousTaxon]
  uniqueTaxa <- unique(taxa)

  tokens <- sub(taxonLine.pattern, "\\5", matrixLines, perl = TRUE)
  if (continuous) {
    tokens <- strsplit(tokens, "\\s+")
    lengths <- vapply(tokens, length, 0L)
    if (length(unique(lengths)) != 1) {
      stop("Different numbers of tokens in different taxa: ",
           paste(lengths, collapse = ', '))
    }
    tokens <- t(vapply(tokens, I, tokens[[1]]))
  } else {
    tokens <- gsub("\t", "", gsub(" ", "", tokens, fixed = TRUE), fixed = TRUE)
    tokens <- vapply(uniqueTaxa,
                     function (taxon) paste0(tokens[taxa == taxon], collapse = ''),
                     character(1))
    tokens <- NexusTokens(tokens, character_num = character_num, session = session)
  }

  rownames(tokens) <- uniqueTaxa

  # Return:
  tokens
}

#' @param tokens Vector of character strings corresponding to phylogenetic
#'  tokens.
#' @describeIn ExtractTaxa Converts a Nexus string to a vector of character
#'  states.
#' @return `NexusTokens` returns a character vector in which each entry
#' corresponds to the states of a phylogenetic character.
#' @examples
#' NexusTokens('01[01]-?')
#' @export
NexusTokens <- function (tokens, character_num=NULL, session=NULL) {
  tokens.pattern <- "\\([^\\)]+\\)|\\[[^\\]]+\\]|\\{[^\\}]+\\}|\\S"
  matches <- gregexpr(tokens.pattern, tokens, perl=TRUE)

  nChar <- length(matches[[1]])

  if (!is.null(session) && requireNamespace('shiny', quietly = TRUE)) {
    shiny::updateNumericInput(session, 'character_num', max = nChar)
  }

  if (!exists("character_num") || is.null(character_num)) {
    character_num <- seq_len(nChar)
  } else if (any(character_num > nChar) || any(character_num < 1)) {
    return(list("Character number must be between 1 and ", nChar, "."))
    character_num[character_num < 1] <- 1
    character_num[character_num > nChar] <- nChar
  }

  tokens <- t(vapply(regmatches(tokens, matches),
                     function (x) x[character_num, drop=FALSE],
                     character(length(character_num))))
  if (length(character_num) == 1) {
    tokens <- t(tokens)
  } else if (length(character_num) == 0) {
    stop("No characters selected")
  }

  # Return:
  tokens
}

#' Read characters from Nexus file
#'
#' Parses a Nexus file, reading character states and names.
#'
#' Tested with nexus files downloaded from MorphoBank with the "no notes"
#' option, but should also work more generally.
#'
#' Please [report](https://github.com/ms609/TreeTools/issues/new?title=Error+parsing+Nexus+file&body=<!--Tell+me+more+and+attach+your+file...-->)
#' incorrectly parsed files.
#'
#' @param filepath character string specifying location of file
#' @template characterNumParam
#' @template sessionParam
#'
#' @return A matrix whose row names correspond to tip labels, and column names
#'         correspond to character labels, with the attribute `state.labels`
#'         listing the state labels for each character; or a character string
#'         explaining why the character cannot be returned.
#'
#' @template MRS
#' @references
#'   \insertRef{Maddison1997}{TreeTools}
#'
#' @examples
#' fileName <- paste0(system.file(package='TreeTools'),
#'                    '/extdata/input/dataset.nex')
#' ReadCharacters(fileName)
#'
#' @export
ReadCharacters <- function (filepath, character_num = NULL, session = NULL) {

  lines <- readLines(filepath, warn = FALSE) # Missing EOL is quite common, so
                                             # warning not helpful
  nexusComment.pattern <- "\\[[^\\]*\\]"
  lines <- gsub(nexusComment.pattern, "", lines)
  lines <- trimws(lines)
  lines <- lines[lines != ""]

  semicolons <- which(RightmostCharacter(lines) == ';')
  upperLines <- toupper(lines)

  continuous <- length(grep('DATATYPE[\\S\\=]+CONTINUOUS', upperLines)) > 0
  matrixStart <- which(upperLines == 'MATRIX')
  if (length(matrixStart) == 0) {
    return(list("MATRIX block not found in Nexus file."))
  } else if (length (matrixStart) > 1) {
    return(list("Multiple MATRIX blocks found in Nexus file."))
  } else {
    matrixEnd <- semicolons[semicolons > matrixStart][1]
    if (lines[matrixEnd] == ';') matrixEnd <- matrixEnd - 1

    matrixLines <- lines[(matrixStart + 1):matrixEnd]
    tokens <- ExtractTaxa(matrixLines, character_num = character_num,
                          session = session, continuous = continuous)
    if (is.null(character_num)) character_num <- seq_len(ncol(tokens))

    ## Written with MorphoBank format in mind: each label on separate line,
    ## each character introduced by integer and terminated with comma.
    labelStart <- which(upperLines == 'CHARLABELS')
    if (length(labelStart) == 1) {
      labelEnd <- semicolons[semicolons > labelStart][1]
      if (lines[labelEnd] == ';') labelEnd <- labelEnd - 1
      #attr(dat, 'char.labels')
      colnames(tokens) <- lines[labelStart + character_num]
    } else {
      if (length(labelStart) > 1)
        return(list("Multiple CharLabels blocks in Nexus file."))
    }

    stateStart <- which(upperLines == 'STATELABELS')
    if (length(stateStart) == 1) {
      stateEnd <- semicolons[semicolons > stateStart][1]
      stateLines <- lines[stateStart:stateEnd]
      stateStarts <- grep("^\\d+", stateLines)
      stateEnds <- grep("[,;]$", stateLines)
      if (length(stateStarts) != length(stateEnds)) {
        warning("Could not parse character states.")
      } else {
        attr(tokens, 'state.labels') <-
          lapply(character_num, function (i)
            stateLines[(stateStarts[i] + 1):(stateEnds[i] - 1)]
          )
      }
    } else {
      if (length(labelStart) > 1) {
        return(list("Multiple StateLabels blocks in Nexus file."))
      }
    }
  }

  # Return:
  tokens
}


#' @describeIn ReadCharacters Read characters from TNT file
#' @export
ReadTntCharacters <- function (filepath, character_num=NULL, session=NULL) {

  lines <- readLines(filepath, warn=FALSE) # Missing EOL might occur in user-
                                           # generated file, so warning not helpful
  tntComment.pattern <- "'[^']*']"
  lines <- gsub(tntComment.pattern, "", lines)
  lines <- trimws(lines)
  lines <- lines[lines != ""]

  semicolons <- which(RightmostCharacter(lines) == ';')
  upperLines <- toupper(lines)

  matrixStart <- which(upperLines == '&[NUM]')
  if (length(matrixStart) == 0) {
    return(list("&[num] entry not found in TNT file."))
  } else if (length (matrixStart) > 1) {
    return(list("Multiple &[num] entries found in TNT file."))
  } else {
    matrixEnd <- semicolons[semicolons > matrixStart][1]
    if (lines[matrixEnd] == ';') matrixEnd <- matrixEnd - 1

    matrixLines <- lines[(matrixStart + 1):matrixEnd]
    tokens <- ExtractTaxa(matrixLines, character_num, session)
    labelStart <- which(upperLines == 'CHARLABELS')
    if (length(labelStart) == 1) {
      labelEnd <- semicolons[semicolons > labelStart][1]
      if (lines[labelEnd] == ';') labelEnd <- labelEnd - 1
      #attr(dat, 'char.labels')
      colnames(tokens) <- lines[labelStart + character_num]
    } else {
      if (length(labelStart) > 1)
        return(list("Multiple CharLabels blocks in Nexus file."))
    }


    labelStart <- which(upperLines == 'CNAMES')
    if (length(labelStart) == 1) {
      labelEnd <- semicolons[semicolons > labelStart][1]
      if (lines[labelEnd] == ';') labelEnd <- labelEnd - 1
      charLines <- lines[labelStart + character_num]
      charLine.pattern <- "^\\S+\\s\\d+\\s(\\w+)(.*)\\s*;\\s*$"

      # Character labels
      charNames <- gsub(charLine.pattern, "\\1", charLines, perl=TRUE)
      colnames(tokens) <- gsub("_", " ", charNames, fixed=TRUE)

      # State labels
      stateNames <- gsub(charLine.pattern, "\\2", charLines, perl=TRUE)
      attr(tokens, 'state.labels') <- lapply(stateNames, function (line) {
        states <- strsplit(trimws(line), "\\s+", perl=TRUE)[[1]]
        trimws(gsub("_", " ", states, fixed=TRUE))
      })
    } else {
      if (length(labelStart) > 1)
        return(list("Multiple cnames entries in TNT file."))
    }
  }

  # Return:
  tokens
}

#' Matrix to `phyDat`
#'
#' Converts a matrix of tokens to a `phyDat` object.
#'
#' @param tokens matrix of tokens, probably created with [`ReadCharacters`]
#'               or [`ReadTntCharacters`]. Row names should correspond to tip
#'               labels; column names may optionally correspond to
#'               character labels.
#'
#' @return an object of class `phyDat`.
#'
#' @template MRS
#' @keywords internal
#' @export
MatrixToPhyDat <- function (tokens) {
  allTokens <- unique(as.character(tokens))
  tokenNumbers <- seq_along(allTokens)
  names(tokenNumbers) <- allTokens
  matches <- gregexpr("[\\d\\-\\w]", allTokens, perl=TRUE)
  whichTokens <- regmatches(allTokens, matches)
  levels <- sort(unique(unlist(whichTokens)))
  whichTokens[allTokens == '?'] <- list(levels)
  contrast <- 1 * t(vapply(whichTokens, function (x) levels %in% x,
                           logical(length(levels))))
  rownames(contrast) <- allTokens
  colnames(contrast) <- levels
  dat <- phangorn::phyDat(tokens, type='USER', contrast=contrast)

  # Return:
  dat
}

#' @describeIn MatrixToPhyDat Converts a `phyDat` object to a matrix of tokens.
#' @param dataset A dataset of class `phyDat`.
## @param parentheses Character vector specifying style of parentheses
## with which to enclose ambiguous characters, e.g, `c('[', ']')` will render
## `[01]`.
## @param sep Character with which to separate ambiguous tokens, e.g. `','`
## will render `[0,1]`.
#' @return `PhyDatToMatrix()` returns a matrix corresponding to the
#' uncompressed character states within a `phyDat` object.
#' @export
PhyDatToMatrix <- function (dataset) {#}, parentheses = c('[', ']'), sep = '') {
  at <- attributes(dataset)
  index <- at$index
  allLevels <- at$allLevels
  t(vapply(dataset, function (x) allLevels[x[index]], character(length(index))))
}

#' @describeIn ReadCharacters Read Nexus characters as `phyDat` object.
#' @importFrom phangorn phyDat
#' @export
ReadAsPhyDat <- function (filepath) {
  MatrixToPhyDat(ReadCharacters(filepath))
}


#' @describeIn ReadCharacters Read TNT characters as `phyDat` object.
#' @importFrom phangorn phyDat
#' @export
ReadTntAsPhyDat <- function (filepath) {
  MatrixToPhyDat(ReadTntCharacters(filepath))
}


#' @describeIn ReadCharacters A convenient wrapper for \pkg{phangorn}'s
#' `phyDat()`, which converts a *list* of morphological characters into a
#' `phyDat` object.
#' If your morphological characters are in the form of a *matrix*, perhaps
#' because they have been read using `read.table`, try [`MatrixToPhyDat()`]
#' instead.
#'
#' @param dataset list of taxa and characters, in the format produced by [read.nexus.data]:
#'                a list of sequences each made of a single vector of mode character,
#'                and named with the taxon name.
#'
#' @export
PhyDat <- function (dataset) {
  nChar <- length(dataset[[1]])
  if (nChar == 1) {
    mat <- matrix(unlist(dataset), dimnames=list(names(dataset), NULL))
  } else {
    mat <- t(vapply(dataset, I, dataset[[1]]))
  }
  MatrixToPhyDat(mat)
}

#' String to `phyDat`
#'
#' Converts a character string to a `phyDat` object.
#'
#' @param string String of tokens, optionally containing whitespace, with no
#'   terminating semi-colon.
#' @param tips Character vector corresponding to the names (in order)
#' of each taxon in the matrix.
#' @param byTaxon Logical; if `TRUE`, string is one **taxon's** coding at a
#' time; if `FALSE`, string is interpreted as one ***character's** coding at a
#' time.
#'
#' @examples
#' morphy <- StringToPhyDat("-?01231230?-", c('Lion', 'Gazelle'), byTaxon = TRUE)
#' # encodes the following matrix:
#' # Lion     -?0123
#' # Gazelle  1230?-
#'
#' @template returnPhydat
#' @seealso \code{\link{phyDat}}
#'
#' @template MRS
#' @aliases StringToPhydat
#' @seealso
#' - [`PhyToString`]
#' @importFrom phangorn phyDat
#' @export
StringToPhyDat <-
  function (string, tips, byTaxon = TRUE) {
    tokens <- matrix(NexusTokens(string), nrow = length(tips), byrow = byTaxon)
    rownames(tokens) <- tips
    MatrixToPhyDat(tokens)
  }
#' @rdname StringToPhyDat
StringToPhydat <- StringToPhyDat

#' `phyDat` to String
#'
#' Extract character data from a `phyDat` object as a string.
#'
#' @param phy An object of class \code{\link{phyDat}}
#' @param parentheses Character specifying format of parentheses with which
#' to surround ambiguous tokens.  Choose from: `\{` (default), `[`, `(`, `<`.
#' @param collapse Character specifying text, perhaps `,`, with which to
#' separate multiple tokens within parentheses
#' @param ps Character specifying text, perhaps `;`, to append to the end of the string
#' @param useIndex Logical (default: `TRUE`) specifying whether to print
#' duplicate characters multiple times, as they appeared in the original matrix
#' @param byTaxon Logical. If `TRUE`, write one taxon followed by the next.
#'                If `FALSE`, write one character followed by the next.
#' @param concatenate Logical specifying whether to concatenate all characters/taxa
#'                    into a single string, or to return a separate string
#'                    for each entry.
#'
#' @examples
#' fileName <- paste0(system.file(package='TreeTools'),
#'                    '/extdata/input/dataset.nex')
#' phyDat <- ReadAsPhyDat(fileName)
#' PhyToString(phyDat, concatenate = FALSE)
#'
#' @return `PhyToString` returns a character vector listing a text
#' representation of the phylogenetic character state for each taxon in turn.
#'
#' @seealso
#' - [`StringToPhyDat`]
#' @template MRS
#' @importFrom phangorn phyDat
#' @export
PhyToString <- function (phy, parentheses = '{', collapse = '', ps = '',
                         useIndex = TRUE, byTaxon = TRUE, concatenate = TRUE) {
  at <- attributes(phy)
  phyLevels <- at$allLevels
  if (sum(phyLevels == '-') > 1) {
    stop("More than one inapplicable level identified.  Is phy$levels malformed?")
  }
  phyChars <- at$nr
  phyContrast <- at$contrast == 1
  phyIndex <- if (useIndex) at$index else seq_len(phyChars)
  outLevels <- at$levels
  inappLevel <- outLevels == '-'

  levelLengths <- vapply(outLevels, nchar, integer(1))
  longLevels <- levelLengths > 1
  if (any(longLevels)) {
    if ('10' %in% outLevels && !(0 %in% outLevels)) {
      outLevels[outLevels == '10'] <- '0'
      longLevels['10'] <- FALSE
    }
    outLevels[longLevels] <- LETTERS[seq_len(sum(longLevels))]
  }

  switch(parentheses,
         '(' = {openBracket <- '('; closeBracket = ')'},
         ')' = {openBracket <- '('; closeBracket = ')'},
         '<' = {openBracket <- '<'; closeBracket = '>'},
         '>' = {openBracket <- '<'; closeBracket = '>'},
         '[' = {openBracket <- '['; closeBracket = ']'},
         ']' = {openBracket <- '['; closeBracket = ']'},
         {openBracket <- '{'; closeBracket = '}'})

  levelTranslation <- apply(phyContrast, 1, function (x)
    ifelse(sum(x) == 1, as.character(outLevels[x]),
           paste0(c(openBracket, paste0(outLevels[x], collapse=collapse),
                    closeBracket), collapse=''))
  )
  if (any(ambigToken <- apply(phyContrast, 1, all))) {
    levelTranslation[ambigToken] <- '?'
  }
  ret <- vapply(phy,
                function (x) levelTranslation[x[phyIndex]],
                character(length(phyIndex)))
  ret <- if (concatenate || is.null(dim(ret))) { # If only one row, don't need to apply
    if (!byTaxon) ret <- t(ret)
    paste0(c(ret, ps), collapse='')
  } else {
    if (byTaxon) ret <- t(ret)
    paste0(apply(ret, 1, paste0, collapse=''), ps)
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
#' @param string Character string.
#' @param len Integer specifying number of characters in `string`.
#'
#' @return Rightmost character of `string`.
#' @examples
#' RightmostCharacter("Hello, World!")
#'
#' @template MRS
#' @export
#' @keywords internal
RightmostCharacter <- function (string, len=nchar(string)) {
  substr(string, len, len)
}

#' Write Newick Tree
#'
#' Writes a tree in Newick format.  This differs from ape's `write.tree`
#' in the encoding of spaces as spaces, rather than underscores.
#'
#' @template treeParam
#'
#' @return `NewickTree` returns a character string denoting `tree` in Newick
#' format.
#'
#' @examples
#' NewickTree(BalancedTree(6))
#'
#' @seealso [`as.Newick`]
#' @importFrom ape write.tree
#' @export
NewickTree <- function(tree) gsub('_', ' ', write.tree(tree), fixed = TRUE)
