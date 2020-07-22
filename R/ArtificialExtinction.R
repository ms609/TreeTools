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
#' @param template Taxon to use as a template.
#' @param subject Vector identifying subject taxa, by name or index.
#' @param replacement Character specifying whether deleted tokens should be
#' replaced with:
#'  - `?`: The ambiguous token, `?`;
#'  - `binary`: The tokens `0` or `1`, with equal probability;
#'  - `uniform`: One of the tokens present in `sampleFrom`, with equal
#'  probability;
#'  - `sample`: One of the tokens present in `sampleFrom`, sampled according
#'  to their frequency.
#' @return A dataset with the same class as `dataset` in which entries that
#' are ambiguous in `template` are made ambiguous in `subject`.
#' @examples
#' set.seed(1)
#' dataset <- matrix(c(sample(0:2, 4 * 8, TRUE),
#'                     '0', '0', rep('?', 6)), nrow = 5,
#'                     dimnames = list(c(LETTERS[1:4], 'FOSSIL'),
#'                                     paste('char', 1:8)), byrow = TRUE)
#' artex <- ArtificialExtinction(dataset, 'FOSSIL', c('A', 'C'))
#' @template MRS
#' @export
ArtificialExtinction <- function (dataset, template, subject,
                                  replacement = '?', sampleFrom = NULL) {
  UseMethod('ArtificialExtinction')
}

#' @rdname ArtificialExtinction
#' @export
ArtificialExtinction.matrix <- function (dataset, template, subject,
                                         replacement = '?',
                                         sampleFrom = NULL) {
  replacers <- c('?', 'binary', 'uniform', 'frequency')
  replaceMode <- pmatch(replacement, replacers)
  if (is.na(replaceMode)) stop("`replacement` unambiguously matched.")

  removes <- dataset[template, ] == '?'
  nRemoves <- sum(removes)
  if (is.null(sampleFrom)) {
    sampleFrom <-
    if (is.numeric(template)) {
      -template
    } else {
      rownames(dataset)[!rownames(dataset) %in% template]
    }
  }

  replaceWith <- switch(replaceMode,
    '?',
    sample(c('0', '1'), nRemoves * length(subject), replace = TRUE), # binary
    apply(unique(dataset[sampleFrom, removes, drop = FALSE]), 2, sample,
          length(subject), replace = TRUE), # Uniform
    apply(dataset[sampleFrom, removes, drop = FALSE], 2, sample,
          length(subject), replace = TRUE) # Frequency
  )
  dataset[subject, removes] <- replaceWith
  dataset
}

#' @rdname ArtificialExtinction
#' @export
ArtificialExtinction.phyDat <- function (dataset, template, subject,
                                         replacement = '?', sampleFrom = NULL) {
  MatrixToPhyDat(ArtificialExtinction(PhyDatToMatrix(dataset), template,
                                      subject, replacement))
}

#' @rdname ArtificialExtinction
#' @export
ArtEx <- ArtificialExtinction
