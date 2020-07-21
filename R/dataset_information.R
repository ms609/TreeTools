#' Mutual information between characters
#' @references \insertRef{Dunn2008}{TreeTools}
#' @examples
#' data('Lobo', package = 'TreeTools')
#' dat <- PhyDatToMatrix(Lobo.phy)[, 1:6]
#'
#' CharacterEntropies(dat)
#' apply(dat, 2, CharacterEntropy)
#' JointCharacterEntropies(dat)
#' MutualCharacterInformation(dat)
#'
#'
#' @template MRS
#' @name CharacterMI
#' @export




#' @rdname CharacterMI
#' @param char,char1,char2 Character vector listing character tokens.
#' For handling of ambiguous tokens, see section 'Ambiguity'.
#' @param ignore Character vector listing tokens to treat as ambiguous,
#' besides `?` and `-`.
#' @section Ambiguity:
#' Ambiguous characters should be specified as a string containing all possible
#' characters states (and, optionally, parentheses): e.g. `[01]`.
#' `?` and `-` will be interpreted as 'any possible value'.
#' Each token is assigned a probability according to its number of
#' *non-ambiguous* occurrences in the character.
#' Example: in `000000 1111 [01] [02] ?`, `[01]` and `?` will be taken to have
#' a 6/10 chance of being token `0`;
#' `[02]` will be treated as having a 100% chance of being state `0`.
#' @examples
#' char <- c(rep('0', 6), rep('1', 4), '[01]', '[02]', '?')
#' ShareAmbiguity(char)
#' ShareAmbiguity(char, ignore = '[02]')
#' CharacterEntropy(char)
#' @export
CharacterEntropy <- function (char, ignore = character(0)) {
  FrequencyToEntropy(ApportionAmbiguity(char))
}

#' @rdname CharacterMI
#' @export
CharacterEntropies <- function (chars, ignore = character(0)) {
  apply(chars, 2, CharacterEntropy, ignore = ignore)
}

#' @rdname CharacterMI
#' @export
ApportionAmbiguity <- function (char, ignore = character(0)) {
  char <- char[!char %in% c('?', '-', ignore)]
  tab <- table(char)
  tokens <- names(table(char))
  singles <- nchar(tokens) == 1L
  unambiguous <- tokens[singles]
  nUnambiguous <- sum(tab[unambiguous])

  ambigs <- vapply(tokens[!singles],
                   function (string) {
                     ret <- tab[unambiguous] * 0
                     x <- tab[strsplit(string, '')[[1]]]
                     x <- tab[string] * x / sum(x, na.rm = TRUE)
                     ret[unambiguous] <- x[unambiguous]
                     ret
                   },
                   tab[unambiguous] * 1)

  ambigSum <- if (!is.null(dim(ambigs))) rowSums(ambigs, na.rm = TRUE) else 0L

  # Return:
  tab[unambiguous] + ambigSum
}


#' @rdname CharacterMI
#' @export
JointCharacterEntropy <- function (char1, char2, ignore = character(0)) {

  ignore <- c('?', '-', ignore)
  ignore1 <- char1 %in% ignore
  ignore2 <- char2 %in% ignore
  bothIgnored <- ignore1 & ignore2

  char1 <- char1[!bothIgnored]
  char2 <- char2[!bothIgnored]

  ignore1 <- ignore1[!bothIgnored]
  ignore2 <- ignore2[!bothIgnored]

  char1[ignore1] <- '?'
  char2[ignore2] <- '?'

  states1 <- table(char1[!ignore1])
  states2 <- table(char2[!ignore2])
  tokens1 <- sort(unique(char1))
  tokens2 <- sort(unique(char2))

  ambigTokens1 <- tokens1 %in% ignore | nchar(tokens1) != 1L
  ambigTokens2 <- tokens2 %in% ignore | nchar(tokens2) != 1L

  confusion <- vapply(tokens1, function (x1) {
    vapply(tokens2,
           function (x2) sum(char1 == x1 & char2 == x2),
           integer(1))
  }, integer(length(unique(char2))))
  if (is.null(dim(confusion))) {
    confusion <- if (length(tokens1) == 1L) {
      cbind(confusion)
    } else {
      rbind(confusion)
    }
  }

  unambigConfusion <- confusion[!ambigTokens2, !ambigTokens1, drop = FALSE]

  known1 <- colSums(confusion[, !ambigTokens1, drop = FALSE])
  known2 <- rowSums(confusion[!ambigTokens2, , drop = FALSE])

  ambig1 <- vapply(tokens1[ambigTokens1], function (token) {
    couldBe <- if (token == '?') names(known1) else strsplit(token, '')[[1]]
    x <- known1[couldBe[couldBe %in% names(known1)]]
    x <- outer(confusion[!ambigTokens2, token], x / sum(x, na.rm = TRUE))
  }, unambigConfusion * 1)
  if (!is.null(dim(ambig1))) {
    ambig1 <- rowSums(ambig1, dims = 2)
  }

  ambig2 <- vapply(tokens2[ambigTokens2], function (token) {
    couldBe <- if (token == '?') names(known2) else strsplit(token, '')[[1]]
    x <- known2[couldBe[couldBe %in% names(known2)]]
    x <- outer(x / sum(x, na.rm = TRUE), confusion[token, !ambigTokens1])
  }, unambigConfusion * 1)
  if (!is.null(dim(ambig2))) {
    ambig2 <- rowSums(ambig2, dims = 2)
  }

  FrequencyToEntropy(unambigConfusion + ambig1 + ambig2)

}

#' @rdname CharacterMI
#' @export
JointCharacterEntropies <- function (characters, ignore = character(0)) {
  nChar <- ncol(characters)
  n <- matrix(seq_len(nChar), nChar, nChar)

  ret <- mapply(function (i, j) JointCharacterEntropy(characters[, i],
                                                      characters[, j],
                                                      ignore = ignore),
                t(n)[lower.tri(n)], n[lower.tri(n)])

  # Return:
  structure(ret, Size = nChar, Diag = FALSE, Upper = FALSE, class = 'dist')
}


#' @rdname CharacterMI
#' @export
MutualCharacterInformation <- function (characters, ignore = character(0)) {

  .AsDist <- function (x) {
    structure(x[lowerTri], Size = nChar, Diag = FALSE, class = 'dist')
  }

  nChar <- ncol(characters)
  m <- nChar - 1L
  jointDist <- JointCharacterEntropies(characters, ignore = ignore)

  h <- CharacterEntropies(characters)
  meanH <- mean(h)

  joint <- as.matrix(jointDist)
  diag(joint) <- h

  mi <- outer(h, h, '+') - joint
  mi[mi < sqrt(.Machine$double.eps)] <- 0
  lowerTri <- lower.tri(mi)

  charMeanMi <- (colSums(mi) - diag(mi)) / m
  meanMi <- mean(mi[lowerTri])

  meanJoint <- mean(jointDist)
  charJoint <- colSums(joint) - diag(joint) / m

  apc <- outer(charMeanMi, charMeanMi) / meanMi
  mip <- .AsDist(mi - apc)

  asc <- outer(charMeanMi, charMeanMi, '+') - meanMi
  mia <- .AsDist(mi - asc)

  mir <- .AsDist(mi / joint)

  # Return:
  list(MI = mi, APC = .AsDist(apc),
       MIr = mir, Zr = (mir - mean(mir)) / sd(mir),
       MIp = mip, Zp = (mip - mean(mip)) / sd(mip),
       MIa = mia, Za = (mia - mean(mia)) / sd(mia)
       )
}

#' @rdname CharacterMI
#' @export
MutCharInfo <- MutualCharacterInformation

#' Entropy of a table of frequencies
#' @param \dots Frequencies  of observations
#' @return `FrequencyToEntropy()` returns the entropy of the input observations.
#' @template MRS
#' @export
FrequencyToEntropy <- function (...) {
  p <- c(...) / sum (...)
  p <- p[p > 0]
  -sum(p * log2(p))
}
