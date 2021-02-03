#' Artificial Extinction
#'
#' Remove tokens that do not occur in a fossil 'template' taxon from a living
#' taxon, to simulate the process of fossilization in removing data from
#' a phylogenetic dataset.
#'
#' Note: this simple implementation does not account for character contingency,
#' e.g. characters whose absence imposes inapplicable or absent tokens on
#' dependent characters.
#'
#' @param dataset Phylogenetic dataset of class `phyDat` or `matrix`.
#' @param subject Vector identifying subject taxa, by name or index.
#' @param template Character or integer identifying taxon to use as a template.
#' @param replaceAmbiguous,replaceCoded Character specifying whether tokens
#' that are ambiguous (`?`) or coded (not `?`) in the fossil template should
#' be replaced with:
#'  - `original`: Their original value; i.e. no change;
#'  - `ambiguous`: The ambiguous token, `?`;
#'  - `binary`: The tokens `0` or `1`, with equal probability;
#'  - `uniform`: One of the tokens present in `sampleFrom`, with equal
#'  probability;
#'  - `sample`: One of the tokens present in `sampleFrom`, sampled according
#'  to their frequency.
#' @param replaceAll Logical: if `TRUE`, replace all tokens in a subject; if
#' `FALSE`, leave any ambiguous tokens (`?`) ambiguous.
#' @param sampleFrom Vector identifying a subset of characters from which to
#' sample replacement tokens.
#' If `NULL`, replacement tokens will be sampled from the initial states of
#' all taxa not used as a template (including the subjects).
#' @return A dataset with the same class as `dataset` in which entries that
#' are ambiguous in `template` are made ambiguous in `subject`.
#' @examples
#' set.seed(1)
#' dataset <- matrix(c(sample(0:2, 4 * 8, TRUE),
#'                     '0', '0', rep('?', 6)), nrow = 5,
#'                     dimnames = list(c(LETTERS[1:4], 'FOSSIL'),
#'                                     paste('char', 1:8)), byrow = TRUE)
#' artex <- ArtificialExtinction(dataset, c('A', 'C'), 'FOSSIL')
#' @template MRS
#' @export
ArtificialExtinction <- function (dataset, subject, template,
                                  replaceAmbiguous = 'ambig',
                                  replaceCoded = 'original',
                                  replaceAll = TRUE,
                                  sampleFrom = NULL) {
  UseMethod('ArtificialExtinction')
}

#' @rdname ArtificialExtinction
#' @export
ArtificialExtinction.matrix <- function (dataset, subject, template,
                                         replaceAmbiguous = 'ambig',
                                         replaceCoded = 'original',
                                         replaceAll = TRUE,
                                         sampleFrom = NULL) {
  replacers <- c('original', 'ambiguous', 'binary', 'uniform', 'frequency')
  replaceA <- pmatch(replaceAmbiguous, replacers)
  if (is.na(replaceA)) stop("`replaceAmbiguous` ambiguously matched.")
  replaceC <- pmatch(replaceCoded, replacers)
  if (is.na(replaceC)) stop("`replaceCoded` ambiguously matched.")

  removes <- dataset[template, ] == '?'
  if (is.null(sampleFrom)) {
    sampleFrom <-
    if (is.numeric(template)) {
      -template
    } else {
      rownames(dataset)[!rownames(dataset) %in% template]
    }
  }

  .DoReplace <- function (dataset, subject, columns, replace) {
    nCols <- sum(columns)
    replaceWith <- switch(replace,
      dataset[subject, columns], # Original
      '?', # Ambiguous
      sample(c('0', '1'), nCols * length(subject), replace = TRUE), # Binary
      apply(unique(dataset[sampleFrom, columns, drop = FALSE]), 2, sample,
            length(subject), replace = TRUE), # Uniform
      apply(dataset[sampleFrom, columns, drop = FALSE], 2, sample,
            length(subject), replace = TRUE) # Frequency
    )

    # Until require R >= 3.5.0
    isFALSE <- function (x) is.logical(x) && length(x) == 1L && !is.na(x) && !x
    if (isFALSE(replaceAll)) {
      replaceWith[dataset[subject, columns] == '?'] <- '?'
    }
    dataset[subject, columns] <- replaceWith
    dataset
  }

  dataset <- .DoReplace(dataset, subject, removes, replaceA)
  dataset <- .DoReplace(dataset, subject, !removes, replaceC)
  dataset
}

#' @rdname ArtificialExtinction
#' @export
ArtificialExtinction.phyDat <- function (dataset, subject, template,
                                         replaceAmbiguous = 'ambig',
                                         replaceCoded = 'original',
                                         replaceAll = TRUE,
                                         sampleFrom = NULL) {
  MatrixToPhyDat(ArtificialExtinction(PhyDatToMatrix(dataset), subject,
                                      template,
                                      replaceAmbiguous,
                                      replaceCoded, sampleFrom))
}

#' @rdname ArtificialExtinction
#' @export
ArtEx <- ArtificialExtinction
