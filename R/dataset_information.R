#' Mutual information between characters
#'
#' Calculate the entropy, joint entropy and mutual information between
#' characters in a morphological dataset,
#' distinguishing 'background' mutual information due to
#' shared ancestry and random chance from functionally significant correlations.
#'
#' `CharacterEntropy()` and `CharacterEntropies()` calculate the entropy of
#' individual characters, using `ApportionAmiguity()` to assign frequencies
#' to ambiguous tokens.
#'
#' `JointCharacterEntropy()` and `JointCharacterEntropies()` calculate the
#' joint entropy of pairs of characters.
#'
#' `MutualCharacterInformation()` calculates the mutual information between
#' each pair of characters: i.e. the extent to which the state of one character
#' can be guessed based on the state of the other.
#'
#' The mutual information (MI) of two characters is assumed to represent the
#' contribution of two components: background MI (_MIb_), comprising
#' a signal from shared ancestry overprinted with random noise, and
#' MI arising from structure and function (_MSsf_).
#' Dunn _et al._ (2008) use this
#' measure to identify amino acids that play a role in the structure or
#' function of proteins -- mutations in one AA will require compensatory
#' mutations in other AA if the protein is to continue functioning.
#' A morphological analogy is proposed.
#'
#' The average product correlation (APC), the product of the average MI of
#' each character with  each other character in the dataset), approximates
#' the background mutual information.
#' A high mean value of MIb across a dataset implies that there is a strong
#' phylogenetic signal.
#'
#' _MIp_, the difference between the mutual information of two characters
#' and their APC, approximates _MIsf_.
#'
#'
#' The runtime increases with the square of the number of characters --
#' may take a couple of seconds to calculate with 100 characters.
#'
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
#' h1 <- JointCharacterEntropies(dat)
#' fit <- cmdscale(h1, k=2)
#' plot(fit, type = 'n')
#' text(fit[, 1], fit[, 2], 1:196, cex = 0.8)
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
JointCharacterEntropy <- function (char1, char2, ignore = character(0),
                                   ignore1 = NULL, ignore2 = NULL,
                                   states1 = NULL, states2 = NULL,
                                   tokens1 = NULL, tokens2 = NULL) {

  if (!isFALSE(ignore)) {
    ignore <- unique(c('?', '-', ignore))
    ignore1 <- char1 %in% ignore
    ignore2 <- char2 %in% ignore
    bothIgnored <- ignore1 & ignore2
    char1 <- char1[!bothIgnored]
    char2 <- char2[!bothIgnored]

    ignore1 <- ignore1[!bothIgnored]
    ignore2 <- ignore2[!bothIgnored]

    char1[ignore1] <- '?'
    char2[ignore2] <- '?'
  }

  if (is.null(states1)) states1 <- table(char1[!ignore1])
  if (is.null(states2)) states2 <- table(char2[!ignore2])
  if (is.null(tokens1)) tokens1 <- sort(unique(char1))
  if (is.null(tokens2)) tokens2 <- sort(unique(char2))

  ambigTokens1 <- tokens1 %in% ignore | nchar(tokens1) != 1L
  ambigTokens2 <- tokens2 %in% ignore | nchar(tokens2) != 1L

  #confusion <- table(char1, char2) # Doubles runtime cf:
  confusion <- vapply(tokens2, function (x2) {
    vapply(tokens1,
           function (x1) sum(char1 == x1 & char2 == x2),
           integer(1))
  }, integer(length(tokens1)))
  if (is.null(dim(confusion))) {
    confusion <- if (length(tokens1) == 1L) {
      rbind(confusion)
    } else {
      cbind(confusion)
    }
  }

  unambigConfusion <- confusion[!ambigTokens1, !ambigTokens2, drop = FALSE]

  known1 <- rowSums(confusion[!ambigTokens1, , drop = FALSE])
  known2 <- colSums(confusion[, !ambigTokens2, drop = FALSE])

  blankConfusion <- unambigConfusion * 0
  ambig1 <- vapply(tokens1[ambigTokens1], function (token) {
    couldBe <- if (token == '?') names(known1) else strsplit(token, '')[[1]]
    x <- known1[couldBe[couldBe %in% names(known1)]]
    # TODO YOU ARE HERE: What if [23] only occurs with `?`?
    x <- outer(x / sum(x, na.rm = TRUE), confusion[token, !ambigTokens2])
    ret <- blankConfusion
    ret[rownames(x), colnames(x)] <- x
    ret
  }, blankConfusion)
  if (!is.null(dim(ambig1))) {
    ambig1 <- rowSums(ambig1, dims = 2)
  }

  ambig2 <- vapply(tokens2[ambigTokens2], function (token) {
    couldBe <- if (token == '?') names(known2) else strsplit(token, '')[[1]]
    x <- known2[couldBe[couldBe %in% names(known2)]]
    x <- outer(confusion[!ambigTokens1, token], x / sum(x, na.rm = TRUE))
    ret <- blankConfusion
    ret[rownames(x), colnames(x)] <- x
    ret
  }, blankConfusion)
  if (!is.null(dim(ambig2))) {
    ambig2 <- rowSums(ambig2, dims = 2)
  }

  FrequencyToEntropy(unambigConfusion + ambig1 + ambig2)

}


#' @rdname CharacterMI
#' @export
JointCharacterEntropies <- function (characters, ignore = character(0)) {
  ignore <- unique(c('-', '?', ignore))
  nChar <- ncol(characters)
  n <- matrix(seq_len(nChar), nChar, nChar)
  tables <- apply(characters, 2L, table)
  tables <- lapply(tables, function (x) x[!names(x) %in% ignore])
  tokens <- apply(characters, 2L, unique)
  tokens <- lapply(tokens, function (x) unique(ifelse(x %in% ignore, '?', x)))
  tokens <- lapply(tokens, sort)

  ignores <- characters %in% ignore

  ret <- mapply(function (i, j) {
    # message(i, ', ', j)
    JointCharacterEntropy(characters[, i], characters[, j],
                          ignore = FALSE,
                          ignore1 = ignores[, i], ignore2 = ignores[, j],
                          states1 = tables[[i]], states2 = tables[[j]],
                          tokens1 = tokens[[i]], tokens2 = tokens[[j]])
    }, t(n)[lower.tri(n)], n[lower.tri(n)])

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
       MIa = mia, Za = (mia - mean(mia)) / sd(mia),
       H = jointDist
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
