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
#' Ambiguity not yet supported, except `?` / `-`.
#' @param ignore Character vector listing tokens to treat as ambiguous,
#' besides `?` and `-`.
#' @export
CharacterEntropy <- function (char, ignore = character(0)) {
  char <- char[!char %in% c('?', '-', ignore)]
  FrequencyToEntropy(table(char))
}


#' @rdname CharacterNI
#' @param char
JointCharacterEntropy <- function (char1, char2, ignore = character(0)) {
  ignore <- c('?', '-', ignore)
  ignore1 <- char1 %in% ignore
  ignore2 <- char2 %in% ignore
  oneIgnored <- ignore1 | ignore2
  states1 <- unique(char1[!ignore1])
  states2 <- unique(char2[!ignore2])
  tokens1 <- sort(unique(char1))
  tokens2 <- sort(unique(char2))

  confusion <- vapply(tokens1, function (x1) {
    vapply(tokens2,
           function (x2) sum(char1 == x1 & char2 == x2),
           integer(1))
  }, integer(length(unique(char2))))

  confusionKnown <- confusion[tokens1[!tokens1 %in% ignore],
                              tokens2[!tokens2 %in% ignore], drop = FALSE]

  rowSums(confusion[tokens1, ignore[ignore %in% tokens2], drop = FALSE])

  unlist(lapply(states1, function (state) {
    states2 <- table(char2[char1 == state])
    ambiguous <- names(states2) %in% ignore
    set2 <- states2[!ambiguous]
    nSet <- sum(set2)


    set2 + (sum(states2[ambiguous]) * set2 / nSet)
  }))

}

#' Entropy of a table of frequencies
FrequencyToEntropy <- function (...) {
  p <- c(...) / sum (...)
  -sum(p * log2(p))
}
