#' Mutual information between characters
#' @references \insertRef{Dunn2008}{TreeTools}
#' @examples
#' data('Lobo', package = 'TreeTools')
#' dat <- PhyDatToMatrix(Lobo.phy)
#'
#' characterEntropies <- apply(dat, 2, CharacterEntropy)
#' head(characterEntropies)
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

  # Return:
  tab[unambiguous] + rowSums(ambigs, na.rm = TRUE)
}


#' @rdname CharacterMI
JointCharacterEntropy <- function (char1, char2, ignore = character(0)) {

  ignore <- c('?', '-', ignore)
  ignore1 <- char1 %in% ignore
  ignore2 <- char2 %in% ignore
  oneIgnored <- ignore1 | ignore2
  bothIgnored <- ignore1 & ignore2

  char1 <- char1[!bothIgnored]
  char2 <- char2[!bothIgnored]
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

  unambigConfusion <- confusion[!ambigTokens2, !ambigTokens1, drop = FALSE]

  known1 <- colSums(confusion[, !ambigTokens1, drop = FALSE])
  known2 <- rowSums(confusion[!ambigTokens2, , drop = FALSE])

  ambig1 <- vapply(tokens1[ambigTokens1], function (token) {
    couldBe <- if (token == '?') names(known1) else strsplit(token, '')[[1]]
    x <- known1[couldBe[couldBe %in% names(known1)]]
    x <- outer(confusion[!ambigTokens2, token], x / sum(x, na.rm = TRUE))
  }, unambigConfusion * 1)

  ambig2 <- vapply(tokens2[ambigTokens2], function (token) {
    couldBe <- if (token == '?') names(known2) else strsplit(token, '')[[1]]
    x <- known2[couldBe[couldBe %in% names(known2)]]
    x <- outer(x / sum(x, na.rm = TRUE), confusion[token, !ambigTokens1])
  }, unambigConfusion * 1)

  FrequencyToEntropy(unambigConfusion +
                     rowSums(ambig1, dims = 2) +
                     rowSums(ambig2, dims = 2))

}


#' @rdname CharacterMI
CharacterEntropy <- function (characters) {
  nChar <- ncol(characters)
  n <- matrix(seq_len(nChar), nChar, nChar)

  is <- t(n)[lower.tri(n)]
  js <- n[lower.tri(n)]

  char1 <- characters[, 1]
  char2 <- characters[, 8]

  JointCharacterEntropy(char1, char2)

  mapply(function (i, j) JointCharacterEntropy(characters[, i], characters[, j]),
                                               is, js)
}

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
